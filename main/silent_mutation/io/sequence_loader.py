"""
silent_mutation.io.sequence_loader

Build a `Variant` directly from user-provided sequences, WITHOUT any genome
FASTA / transcript lookup / ClinVar VCF. This is the entry point for the
"paste a sequence" web-tool mode.

Two coding-context modes
------------------------
1. CONTIGUOUS (default, no exon args): the whole window is assumed in-frame CDS.
   Frame comes from `var_codon_phase` (or is inferred from the edit position
   when `wt_in_frame_from_start=True`). This is the original demo behaviour.

2. EXON-AWARE (pass exon_start / exon_end / codon_start): the user marks which
   part of the window is exon (coding) and where a codon starts (reading frame).
   codon_lookup is then built ONLY over the exonic codons. This supports the
   key case where the EDIT is intronic (e.g. a splice-acceptor correction
   c.8229-2A>G) but the SILENT mutation must land in the adjacent EXON:
     * intron positions (incl. the edit) are absent from codon_lookup, so
       silent_finder ignores them automatically;
     * the 1-2 nt split-codon "leftover" between the exon edge and the first
       full codon is likewise absent (it needs the neighbouring exon to form a
       codon, which isn't in this local window) and is therefore ignored too.

Two ways to specify the edit
----------------------------
1. build_variant_from_pair(wt, ed)             # paste WT + Edited, we diff them
2. build_variant_from_spec(wt, pos, ref, alt)  # paste WT + an explicit edit

Both funnel into _assemble_variant(), which carves a 60bp-flank window
(VAR_IDX = 60) around the edit, builds seq_ed, builds codon_lookup (contiguous
or exon-aware), and returns a fully-validated Variant.

Coordinate convention
----------------------
genome_pos = VAR_IDX (60) and cds_strand="+", so Variant.seq_idx_to_genome_pos(i)
== i: genome positions and seq_wt indices coincide, and codon_lookup is keyed by
that shared coordinate. The exon_start / exon_end / codon_start arguments are
0-based indices into the WT STRING you pass in (the same coordinates a sequence
viewer / slider bar would produce); they are translated to window coordinates
internally.
"""

from __future__ import annotations

from dataclasses import dataclass
from typing import Optional

from silent_mutation.core.types import VAR_IDX, FLANK, Variant
from silent_mutation.core.codon_utils import load_codon_table


# ─────────────────────────────────────────────────────────────────────────────
# Lightweight CodonInfo (attribute-compatible with io.codon_lookup.CodonInfo)
# ─────────────────────────────────────────────────────────────────────────────

@dataclass(frozen=True)
class CodonInfo:
    codon: str
    frame: int
    aa: str
    codon_index: int
    codon_genomic_positions: tuple[int, int, int]


_VALID = set("ACGT")


def _clean(seq: str, name: str) -> str:
    s = seq.strip().upper().replace(" ", "").replace("\n", "")
    bad = set(s) - _VALID
    if bad:
        raise ValueError(f"{name} contains non-ACGT characters: {sorted(bad)}")
    return s


# ─────────────────────────────────────────────────────────────────────────────
# Edit detection (minimal ref/alt by prefix/suffix trimming)
# ─────────────────────────────────────────────────────────────────────────────

def diff_edit(wt: str, ed: str) -> tuple[int, str, str]:
    """Locate the single contiguous edit between WT and ED by trimming the
    common prefix and common suffix. Returns (edit_pos, ref, alt) where
    edit_pos is the 0-based index in `wt` at which the edit begins.
    Handles 1-3bp substitution / insertion / deletion.
    """
    if wt == ed:
        raise ValueError("WT and Edited sequences are identical — no edit found.")

    n = min(len(wt), len(ed))
    i = 0
    while i < n and wt[i] == ed[i]:
        i += 1
    j = 0
    while j < n - i and wt[len(wt) - 1 - j] == ed[len(ed) - 1 - j]:
        j += 1

    ref = wt[i:len(wt) - j]
    alt = ed[i:len(ed) - j]

    if len(ref) > 3 or len(alt) > 3:
        raise ValueError(
            f"Edit is {len(ref)}bp→{len(alt)}bp; only 1-3bp edits are supported. "
            f"(ref={ref!r}, alt={alt!r})"
        )
    if len(ref) == 0 and len(alt) == 0:
        raise ValueError("Empty edit after trimming — sequences may differ oddly.")
    return i, ref, alt


# ─────────────────────────────────────────────────────────────────────────────
# Contiguous codon_lookup (no introns — in-frame assumption)  [UNCHANGED]
# ─────────────────────────────────────────────────────────────────────────────

def _build_codon_lookup(
    seq_wt: str,
    var_codon_phase: int,
    var_codon_number: Optional[int],
    table: dict[str, str],
) -> tuple[dict[int, CodonInfo], int, int]:
    """Build a contiguous codon_lookup over the whole window (in-frame CDS)."""
    L = len(seq_wt)
    var_codon_start = VAR_IDX - var_codon_phase          # codon boundary holding VAR_IDX
    first_start = var_codon_start % 3                     # first codon start in window
    starts = list(range(first_start, L - 2, 3))
    anchor = starts.index(var_codon_start) if var_codon_start in starts else None

    lookup: dict[int, CodonInfo] = {}
    for k, s in enumerate(starts):
        codon = seq_wt[s:s + 3]
        if var_codon_number is not None and anchor is not None:
            cidx = var_codon_number + (k - anchor)
        else:
            cidx = k + 1
        aa = table.get(codon, "?")
        gpos = (s, s + 1, s + 2)                          # genome_pos == seq_idx
        for fr in range(3):
            lookup[s + fr] = CodonInfo(codon, fr, aa, cidx, gpos)
    return lookup, first_start, len(starts)


# ─────────────────────────────────────────────────────────────────────────────
# Exon-aware codon_lookup (single exonic segment of the window)  [NEW]
# ─────────────────────────────────────────────────────────────────────────────

def _build_codon_lookup_exon(
    seq_wt: str,
    exon_start_w: int,
    exon_end_w: int,
    codon_start_w: int,
    base_codon_no: Optional[int],
    table: dict[str, str],
) -> tuple[dict[int, CodonInfo], int, int]:
    """Build codon_lookup over ONE exonic segment of the window only.

    All indices are WINDOW (seq_wt) coordinates.
        coding region : [exon_start_w, exon_end_w)
        reading frame : anchored on codon_start_w  (only codon_start_w % 3 and
                        the codon numbering depend on it)
        base_codon_no : codon number of the codon at codon_start_w (true protein
                        numbering); None -> number locally from 1.

    Codons start at the first in-frame boundary >= exon_start_w and step by 3,
    each fully inside [exon_start_w, exon_end_w). Everything else — intron
    (outside the exon) and the split-codon leftover between the exon edge and the
    first full codon — is simply absent from the lookup, so silent_finder ignores
    it automatically. The edit at VAR_IDX, if intronic, is naturally not present.
    """
    L = len(seq_wt)
    e_lo = max(exon_start_w, 0)
    e_hi = min(exon_end_w, L)
    frame = codon_start_w % 3                              # codon boundaries: p % 3 == frame
    first = e_lo + ((frame - e_lo) % 3)                    # first in-frame boundary >= e_lo

    lookup: dict[int, CodonInfo] = {}
    starts = list(range(first, e_hi - 2, 3))               # each codon fully within the exon
    for k, s in enumerate(starts):
        codon = seq_wt[s:s + 3]
        if base_codon_no is not None:
            cidx = base_codon_no + (s - codon_start_w) // 3   # both ≡ frame (mod 3) → exact
        else:
            cidx = k + 1
        aa = table.get(codon, "?")
        gpos = (s, s + 1, s + 2)
        for fr in range(3):
            lookup[s + fr] = CodonInfo(codon, fr, aa, cidx, gpos)
    first_start = starts[0] if starts else first
    return lookup, first_start, len(starts)


# ─────────────────────────────────────────────────────────────────────────────
# Core assembler
# ─────────────────────────────────────────────────────────────────────────────

def _assemble_variant(
    wt: str,
    edit_pos: int,
    ref: str,
    alt: str,
    var_codon_phase: Optional[int],
    var_codon_number: Optional[int] = None,
    variant_id: str = "user_seq",
    gene_symbol: str = "",
    table: Optional[dict[str, str]] = None,
    *,
    exon_start: Optional[int] = None,
    exon_end: Optional[int] = None,
    codon_start: Optional[int] = None,
    cds_coord_at_codon_start: Optional[int] = None,
) -> Variant:
    """Carve the 60bp-flank window and build a fully-validated Variant.

    Exon args (exon_start / exon_end / codon_start) are in WT-input coordinates
    and are translated to window coordinates here. If any one is supplied, all
    three are required.
    """
    table = table or load_codon_table()
    ref_len, alt_len = len(ref), len(alt)

    lo = edit_pos - FLANK
    hi = edit_pos + ref_len + FLANK
    if lo < 0 or hi > len(wt):
        raise ValueError(
            f"Not enough flanking sequence: need >= {FLANK}bp on each side of the "
            f"edit (edit at index {edit_pos} in a {len(wt)}bp WT). "
            f"Provide a longer WT sequence."
        )

    seq_wt = wt[lo:hi]                                    # FLANK + ref_len + FLANK
    seq_ed = seq_wt[:VAR_IDX] + alt + seq_wt[VAR_IDX + ref_len:]

    exon_mode = (exon_start is not None) or (exon_end is not None) or (codon_start is not None)
    if exon_mode:
        if exon_start is None or exon_end is None or codon_start is None:
            raise ValueError(
                "Exon mode needs exon_start, exon_end and codon_start together "
                "(0-based indices into the WT string)."
            )
        if not (exon_start <= codon_start and exon_start < exon_end):
            raise ValueError(
                f"Expect exon_start <= codon_start and exon_start < exon_end "
                f"(got exon_start={exon_start}, exon_end={exon_end}, "
                f"codon_start={codon_start})."
            )
        # WT-input coords -> window coords
        es_w = exon_start - lo
        ee_w = exon_end - lo
        cs_w = codon_start - lo
        base_codon_no = (
            (cds_coord_at_codon_start - 1) // 3 + 1
            if cds_coord_at_codon_start is not None else None
        )
        lookup, first_start, n_codons = _build_codon_lookup_exon(
            seq_wt, es_w, ee_w, cs_w, base_codon_no, table,
        )
    else:
        lookup, first_start, n_codons = _build_codon_lookup(
            seq_wt, var_codon_phase, var_codon_number, table,
        )

    return Variant(
        seq_wt=seq_wt,
        seq_ed=seq_ed,
        ref_len=ref_len,
        alt_len=alt_len,
        cds_strand="+",
        cds_frame=0,                                     # cds_start is a codon boundary
        cds_start_in_seq=first_start,
        cds_end_in_seq=first_start + n_codons * 3,
        variant_id=variant_id,
        gene_symbol=gene_symbol,
        transcript_id="",
        chrom="",
        genome_pos=VAR_IDX,                              # makes genome_pos == seq_idx
        ref_allele=ref,                                  # validated in __post_init__
        alt_allele=alt,
        codon_lookup=lookup,
    )


# ─────────────────────────────────────────────────────────────────────────────
# Public entry points
# ─────────────────────────────────────────────────────────────────────────────

def build_variant_from_pair(
    wt: str,
    ed: str,
    var_codon_phase: Optional[int] = None,
    wt_in_frame_from_start: bool = True,
    var_codon_number: Optional[int] = None,
    variant_id: str = "user_seq",
    gene_symbol: str = "",
    exon_start: Optional[int] = None,
    exon_end: Optional[int] = None,
    codon_start: Optional[int] = None,
    cds_coord_at_codon_start: Optional[int] = None,
) -> Variant:
    """Build a Variant from a WT + Edited sequence pair (we diff them).

    Contiguous mode (no exon args): frame from var_codon_phase, or inferred from
    the edit position when wt_in_frame_from_start=True.

    Exon-aware mode (pass exon_start / exon_end / codon_start, 0-based indices
    into `wt`): codon_lookup is built only over the exon, with the reading frame
    anchored on codon_start. Use this when the edit is intronic but silent
    markers must land in the adjacent exon. var_codon_phase is ignored here —
    the frame comes from codon_start. Optionally give cds_coord_at_codon_start
    (the c. coordinate of the codon_start base) for true protein numbering.
    """
    wt = _clean(wt, "WT")
    ed = _clean(ed, "Edited")
    edit_pos, ref, alt = diff_edit(wt, ed)

    exon_mode = (exon_start is not None) or (exon_end is not None) or (codon_start is not None)
    if not exon_mode and var_codon_phase is None:
        if not wt_in_frame_from_start:
            raise ValueError(
                "var_codon_phase is required when WT does not start at a codon boundary."
            )
        var_codon_phase = edit_pos % 3

    return _assemble_variant(
        wt, edit_pos, ref, alt, var_codon_phase,
        var_codon_number=var_codon_number,
        variant_id=variant_id, gene_symbol=gene_symbol,
        exon_start=exon_start, exon_end=exon_end, codon_start=codon_start,
        cds_coord_at_codon_start=cds_coord_at_codon_start,
    )


def build_variant_from_spec(
    wt: str,
    pos: int,
    ref: str,
    alt: str,
    var_codon_phase: Optional[int] = None,
    wt_in_frame_from_start: bool = True,
    var_codon_number: Optional[int] = None,
    variant_id: str = "user_seq",
    gene_symbol: str = "",
    exon_start: Optional[int] = None,
    exon_end: Optional[int] = None,
    codon_start: Optional[int] = None,
    cds_coord_at_codon_start: Optional[int] = None,
) -> Variant:
    """Build a Variant from WT + an explicit edit at 0-based index `pos`.
    For substitutions/deletions, wt[pos:pos+len(ref)] must equal ref.
    For insertions, ref="" and the insertion happens before index `pos`.
    Exon args behave exactly as in build_variant_from_pair.
    """
    wt = _clean(wt, "WT")
    ref = _clean(ref, "ref") if ref else ""
    alt = _clean(alt, "alt") if alt else ""

    if ref and wt[pos:pos + len(ref)] != ref:
        raise ValueError(
            f"ref {ref!r} does not match WT at index {pos} "
            f"(found {wt[pos:pos + len(ref)]!r})."
        )
    if len(ref) > 3 or len(alt) > 3:
        raise ValueError("Only 1-3bp edits are supported.")

    exon_mode = (exon_start is not None) or (exon_end is not None) or (codon_start is not None)
    if not exon_mode and var_codon_phase is None:
        if not wt_in_frame_from_start:
            raise ValueError(
                "var_codon_phase is required when WT does not start at a codon boundary."
            )
        var_codon_phase = pos % 3

    return _assemble_variant(
        wt, pos, ref, alt, var_codon_phase,
        var_codon_number=var_codon_number,
        variant_id=variant_id, gene_symbol=gene_symbol,
        exon_start=exon_start, exon_end=exon_end, codon_start=codon_start,
        cds_coord_at_codon_start=cds_coord_at_codon_start,
    )
