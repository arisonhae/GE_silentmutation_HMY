"""
silent_mutation.io.deepprime_runner

DeepPrime efficiency ranking via the genet package (genet.predict.DeepPrime).
Scores candidate pegRNAs and returns them ranked by predicted efficiency.
genet (>= 0.17) reproduces the published SynDesign reference (deepcrispr.info).

run_deepprime_silent(...) returns a list of pegRNA dicts sorted by efficiency;
each dict carries that pegRNA's silent candidates under "outputs". The
/api/analyze use_dp branch in server.py consumes this structure and index.html
renders it. To swap genet for a future co-edit scoring engine, preserve this
input/output contract -- nothing else in server.py / index.html needs to change.

Nick/RTT recovery (build_pam_from_wt74): from genet's WT 74-mer ('Target'
column), the nick sits at Target_pos + 21 on the + strand, reproducing genet's
Edit_pos (e.g. Target@38 -> nick 59 -> edit@60 = Edit_pos 2).

genet is imported lazily. Import errors, and genet's internal sys.exit() on an
unsupported PE-system/cell-type combo, are surfaced as DeepPrimeUnavailable so
the web tool falls back to the standalone candidate list instead of crashing
Flask. Requires Python 3.10 (genet depends on tensorflow<2.10).
"""

from __future__ import annotations

from typing import Optional

import pandas as pd

from silent_mutation.io.sequence_loader import build_variant_from_pair
from silent_mutation.io.genome_loader import reverse_complement
from silent_mutation.core.types import PamCandidate, Variant, VAR_IDX
from silent_mutation.core.silent_finder import find_silent_candidates
from silent_mutation.core.pegrna_builder import build_pegrna_output
from silent_mutation.core.codon_utils import load_codon_table


class DeepPrimeUnavailable(RuntimeError):
    """Raised when genet cannot be imported or run."""


def _edit_notation(variant: Variant) -> str:
    """Build genet's single-sequence input with the edit marked as (ref/alt) at
    VAR_IDX.
        sub: (G/T)   insertion: (/AT)   deletion: (AT/)
    Substitution is verified against genet 0.17. For ins/del the (ref/alt)
    convention follows the natural minimal representation — VERIFY against genet
    before relying on indel scoring.
    """
    i = VAR_IDX
    ref = variant.ref_allele or ""
    alt = variant.alt_allele or ""
    return variant.seq_wt[:i] + f"({ref}/{alt})" + variant.seq_wt[i + variant.ref_len:]


def score_all_pegrnas(
    variant: Variant,
    pe_system: str = "PE2max", cell_type: str = "HEK293T",
    pbs_min: int = 1, pbs_max: int = 17, rtt_max: int = 40,
) -> pd.DataFrame:
    """Run genet DeepPrime and return a df whose columns are renamed to the names
    the rest of this module expects:
        WT74_On (= genet 'Target'), PBSlen, RTlen, DeepPrime_score (= '{pe}_score').
    Sorted by DeepPrime_score desc. Empty df if no pegRNA.
    Raises DeepPrimeUnavailable if genet import/run fails (incl. genet sys.exit()
    on an unsupported pe_system/cell_type combo).
    """
    try:
        from genet.predict import DeepPrime
    except Exception as e:
        raise DeepPrimeUnavailable(
            "genet could not be imported. Run inside the lab 'genet' conda env "
            "(genet 0.17 + torch + tensorflow + ViennaRNA). Original error: %s" % e
        ) from e

    seq = _edit_notation(variant)
    try:
        dp = DeepPrime(seq, name="seq_query", pam="NGG",
                       pbs_min=pbs_min, pbs_max=pbs_max, rtt_min=0, rtt_max=rtt_max,
                       spacer_len=20)
        df = dp.predict(pe_system=pe_system, cell_type=cell_type)
    except (Exception, SystemExit) as e:
        # genet's load_model sys.exit()s on an unsupported PE/cell combo; catch it
        # so an invalid UI selection can't kill the Flask server.
        raise DeepPrimeUnavailable(
            "genet DeepPrime run failed for pe_system=%r cell_type=%r "
            "(this combo may be unsupported by genet): %s"
            % (pe_system, cell_type, e)
        ) from e

    if df is None or len(df) == 0:
        return pd.DataFrame()

    score_col = f"{pe_system}_score"
    if score_col not in df.columns:
        cands = [c for c in df.columns if c.endswith("_score") and "cas9" not in c.lower()]
        if not cands:
            raise DeepPrimeUnavailable(
                "genet output has no score column (cols=%s)" % list(df.columns))
        score_col = cands[0]

    df = df.rename(columns={
        "Target":  "WT74_On",
        "PBS_len": "PBSlen",
        "RTT_len": "RTlen",
        score_col: "DeepPrime_score",
    })

    # genet's 'Edit_pos' is kept as-is. If you want SynDesign-like behaviour
    # (drop nick-proximal candidates that SynDesign filters), uncomment:
    #     df = df[df["Edit_pos"] >= 3]
    return df.sort_values("DeepPrime_score", ascending=False).reset_index(drop=True)


def _spacer_from_pam(seq_wt: str, pam, spacer_len: int = 20) -> str:
    """Protospacer (guide), PAM-strand 5'->3', for a recovered PAM."""
    if pam is None:
        return ""
    if pam.pam_pattern == "NGG":
        s = pam.pam_pos - spacer_len
        return seq_wt[s:pam.pam_pos] if s >= 0 else ""
    s = pam.pam_pos + 3
    return reverse_complement(seq_wt[s:s + spacer_len]) if s + spacer_len <= len(seq_wt) else ""


def build_pam_from_wt74(seq_wt: str, wt74: str, rtlen: int, pbs_len: int) -> Optional[PamCandidate]:
    """Recover nick + strand from the WT 74-mer (genet's 'Target' column) and
    build a PamCandidate whose
    RTT window equals the pegRNA's RTT. Returns None if it can't be placed.

    Geometry verified against genet: nick = Target_pos + 21 (+strand) reproduces
    Edit_pos exactly.
    """
    L = len(seq_wt)

    pos = seq_wt.find(wt74)                # '+' strand: WT74 is a direct substring
    if pos != -1:
        nick = pos + 21
        rtt_start, rtt_end = nick, min(nick + rtlen, L)
        pbs_s = nick - pbs_len
        if pbs_s < 0:
            return None
        return PamCandidate(
            pam_pos=nick + 3, pam_pattern="NGG", nick_pos=nick,
            pbs_seq=seq_wt[pbs_s:nick], rtt_seq=seq_wt[rtt_start:rtt_end],
            rtt_start_in_seq=rtt_start, rtt_end_in_seq=rtt_end,
        )

    pos = seq_wt.find(reverse_complement(wt74))   # '-' strand
    if pos != -1:
        nick = pos + 53
        rtt_end, rtt_start = nick, max(nick - rtlen, 0)
        pbs_e = nick + pbs_len
        if pbs_e > L:
            return None
        return PamCandidate(
            pam_pos=nick - 6, pam_pattern="CCN", nick_pos=nick,
            pbs_seq=reverse_complement(seq_wt[nick:pbs_e]),
            rtt_seq=reverse_complement(seq_wt[rtt_start:rtt_end]),
            rtt_start_in_seq=rtt_start, rtt_end_in_seq=rtt_end,
        )
    return None


def run_deepprime_silent(
    wt: str, ed: str,
    var_codon_phase: Optional[int] = None,
    wt_in_frame_from_start: bool = True,
    var_codon_number: Optional[int] = None,
    pe_system: str = "PE2max", cell_type: str = "HEK293T",
    pbs_min: int = 1, pbs_max: int = 17, rtt_max: int = 40,
    top_n: int = 20, dp_path: str = "",
    exon_start: Optional[int] = None,
    exon_end: Optional[int] = None,
    codon_start: Optional[int] = None,
    cds_coord_at_codon_start: Optional[int] = None,
) -> list:
    """Run genet scoring. Returns an efficiency-ranked list of pegRNA
    dicts; each carries its full list of silent candidates (PegRNAOutput) under
    'outputs'. Return structure matches what server.py and index.html expect.

    `dp_path` is accepted for backward compatibility (server.py still passes it)
    but ignored — genet resolves its own bundled models.

    Exon-aware: pass exon_start / exon_end / codon_start (0-based indices into
    `wt`) and optional cds_coord_at_codon_start to build the codon lookup over the
    exon only. DeepPrime scoring is unaffected (it works on the sequence); the
    silent markers per pegRNA then come from the exon lookup, so an intronic edit
    carries its silent markers in the adjacent exon — same as the standalone path.

    Raises DeepPrimeUnavailable if genet can't run. Returns [] if no pegRNAs.
    """
    variant = build_variant_from_pair(
        wt, ed, var_codon_phase=var_codon_phase,
        wt_in_frame_from_start=wt_in_frame_from_start,
        var_codon_number=var_codon_number,
        exon_start=exon_start, exon_end=exon_end, codon_start=codon_start,
        cds_coord_at_codon_start=cds_coord_at_codon_start,
    )
    df = score_all_pegrnas(
        variant, pe_system=pe_system, cell_type=cell_type,
        pbs_min=pbs_min, pbs_max=pbs_max, rtt_max=rtt_max,
    )
    if df.empty:
        return []

    table = load_codon_table()
    pegrnas: list[dict] = []
    for rank, (_, rec) in enumerate(df.head(top_n).iterrows(), start=1):
        pam = build_pam_from_wt74(
            variant.seq_wt, str(rec["WT74_On"]), int(rec["RTlen"]), int(rec["PBSlen"]),
        )
        outputs = []
        if pam is not None:
            silents = find_silent_candidates(variant, pam, table)
            outputs = [build_pegrna_output(variant, pam, s) for s in silents]
            outputs.sort(key=lambda o: o.priority)        # priority 1 (PAM) first
        best = outputs[0] if outputs else None
        pegrnas.append({
            "rank": rank,
            "efficiency": round(float(rec["DeepPrime_score"]), 2),
            "strand": "+" if (pam and pam.pam_pattern == "NGG") else ("-" if pam else "?"),
            "spacer": _spacer_from_pam(variant.seq_wt, pam),
            "pbs_len": int(rec["PBSlen"]),
            "rtt_len": int(rec["RTlen"]),
            "n_silent": len(outputs),
            "silent_marker": "available" if best else "none",
            "best_codon": f"{best.original_codon}\u2192{best.mutated_codon}" if best else "",
            "best_locale": best.locale if best else "",
            "outputs": outputs,                            # list[PegRNAOutput] for viz
        })
    return pegrnas
