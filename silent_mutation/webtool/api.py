"""
silent_mutation.webtool.api — Flask-free core of the web tool.

This module reproduces the request-handling logic of server.py's /api/analyze
and /api/verify routes as plain functions that take a dict and return a dict,
with NO Flask and NO DeepPrime/genet dependency. It is what the self-contained
(Pyodide) single-file build calls, and server.py can also delegate to it so the
two stay byte-identical.

run_analyze(body) -> dict   (same shape as /api/analyze JSON)
run_verify(body)  -> dict   (same shape as /api/verify JSON)

DeepPrime ranking is intentionally unavailable here (genet is Python-server
only); when body["use_dp"] is set, the response carries
deepprime={"available": False, ...} exactly like the server's graceful
DeepPrimeUnavailable fallback, and index.html already renders that state.
"""

from __future__ import annotations

import traceback

from silent_mutation.io.sequence_loader import (
    build_variant_from_pair, build_variant_from_spec,
)
from silent_mutation.core.pegrna_builder import build_outputs_for_variant
from silent_mutation.io.genome_loader import reverse_complement


# ── visualization: per-base roles (silent > variant > PAM > PBS > RTT > bg) ──
# (verbatim from server.py)

def role_track(row) -> list:
    L = len(row.seq_wt)
    track = ["bg"] * L
    rtt_s, rtt_e = row.rtt_start_in_seq, row.rtt_end_in_seq
    pam_s, pam_e = row.pam_pos, row.pam_pos + 3
    if row.pam_pattern == "NGG":
        pbs_s, pbs_e = row.nick_pos - len(row.pbs_seq), row.nick_pos
    else:
        pbs_s, pbs_e = row.nick_pos, row.nick_pos + len(row.pbs_seq)
    for i in range(max(rtt_s, 0), min(rtt_e, L)):
        track[i] = "rtt"
    for i in range(max(pbs_s, 0), min(pbs_e, L)):
        track[i] = "pbs"
    for i in range(max(pam_s, 0), min(pam_e, L)):
        track[i] = "pam"
    for i in range(60, min(60 + max(row.ref_len, 1), L)):
        track[i] = "variant"
    for p in row.silent_positions:
        if 0 <= p < L:
            track[p] = "silent"
    return track


def _spacer_of(row, spacer_len: int = 20) -> str:
    """Protospacer (guide), PAM-strand 5'->3'."""
    if row.pam_pattern == "NGG":
        s = row.pam_pos - spacer_len
        return row.seq_wt[s:row.pam_pos] if s >= 0 else ""
    s = row.pam_pos + 3
    return reverse_complement(row.seq_wt[s:s + spacer_len]) \
        if s + spacer_len <= len(row.seq_wt) else ""


def _pbsrtt_ext(row):
    """Return (pbsrtt_forward, ext_rc)."""
    try:
        if len(row.seq_ed_with_silent) != len(row.seq_wt):
            return "", ""
        rs, re_ = row.rtt_start_in_seq, row.rtt_end_in_seq
        if row.pam_pattern == "NGG":
            edited_rtt = row.seq_ed_with_silent[rs:re_]
        else:
            edited_rtt = reverse_complement(row.seq_ed_with_silent[rs:re_])
        pbsrtt = row.pbs_seq + edited_rtt
        return pbsrtt, reverse_complement(pbsrtt)
    except Exception:
        return "", ""


def candidate_dict(row, rank: int) -> dict:
    if row.pam_pattern == "NGG":
        sil_off = [p - row.rtt_start_in_seq for p in row.silent_positions]
    else:
        sil_off = [(row.rtt_end_in_seq - 1) - p for p in row.silent_positions]
    pbsrtt, ext_rc = _pbsrtt_ext(row)
    return {
        "rank": rank,
        "pam": f"{row.pam_pattern}@{row.pam_pos}",
        "pam_pattern": row.pam_pattern,
        "locale": row.locale,
        "priority": row.priority,
        "codon": f"{row.original_codon}\u2192{row.mutated_codon}",
        "aa": row.aa,
        "silent_pos": ",".join(map(str, row.silent_positions)),
        "silent": f"{row.silent_ref}\u2192{row.silent_alt}",
        "rtt": row.rtt_seq_wt,
        "spacer": _spacer_of(row),
        "pbs_seq": row.pbs_seq,
        "pbs_len": len(row.pbs_seq),
        "rtt_len": len(row.rtt_seq_wt),
        "pbsrtt": pbsrtt,
        "ext_rc": ext_rc,
        "edseq": row.seq_ed_with_silent,
        "seq_wt": row.seq_wt,
        "roles": role_track(row),
        "nick": row.nick_pos,
        "sil_off": sil_off,
    }


def _exon_kwargs(d) -> dict:
    """Pull exon-aware args from the request body, if the bar sent them."""
    if (d.get("exon_start") is None or d.get("exon_end") is None
            or d.get("codon_start") is None):
        return {}
    cc = str(d.get("cds_coord", "")).strip()
    has_cds = cc.lstrip("-").isdigit()
    return dict(
        exon_start=int(d["exon_start"]),
        exon_end=int(d["exon_end"]),
        codon_start=int(d["codon_start"]),
        cds_coord_at_codon_start=int(cc) if has_cds else None,
    )


# ── analyze ──

def run_analyze(d: dict) -> dict:
    """Same contract as server.py /api/analyze, minus Flask, minus DeepPrime."""
    try:
        phase_map = {"0": 0, "1": 1, "2": 2}
        frame = d.get("frame", "start")
        var_phase = phase_map.get(frame)            # None => infer from boundary
        from_start = frame == "start"
        cn = str(d.get("codon_number", "")).strip()
        var_num = int(cn) if cn.isdigit() else None

        ex = _exon_kwargs(d)

        mode = d.get("mode", "pair")
        if mode == "pair":
            wt, ed = d["wt"], d["ed"]
            variant = build_variant_from_pair(
                wt, ed, var_codon_phase=var_phase,
                wt_in_frame_from_start=from_start, var_codon_number=var_num,
                **ex,
            )
        else:
            wt = d["wt"].strip().upper().replace(" ", "").replace("\n", "")
            pos, ref, alt = int(d["pos"]), d.get("ref", "").upper(), d.get("alt", "").upper()
            variant = build_variant_from_spec(
                wt, pos, ref, alt, var_codon_phase=var_phase,
                wt_in_frame_from_start=from_start, var_codon_number=var_num,
                **ex,
            )
    except Exception as e:
        return {"ok": False, "error": str(e)}

    pbs_len = int(d.get("pbs_len", 13))
    max_rtt = int(d.get("max_rtt", 40))
    rows = build_outputs_for_variant(variant, pbs_len=pbs_len, max_rtt_len=max_rtt)

    info = variant.codon_at_variant()
    protein = rows[0].variant_protein_change if rows else ""
    if not protein and variant.variant_type != "substitution":
        net = variant.alt_len - variant.ref_len
        protein = "frameshift" if (net % 3 != 0) else "in-frame indel"
    resp = {
        "ok": True,
        "variant": {
            "ref": variant.ref_allele or "-",
            "alt": variant.alt_allele or "-",
            "type": variant.variant_type,
            "codon": info.codon if info else "-",
            "aa": info.aa if info else "-",
            "protein": protein,
        },
        "candidates": [candidate_dict(r, i + 1) for i, r in enumerate(rows)],
        "deepprime": None,
    }

    # DeepPrime is server-only (genet). Mirror the server's graceful-unavailable
    # response so index.html renders the standalone list with a clear note.
    if d.get("use_dp"):
        resp["deepprime"] = {
            "available": False,
            "pegrnas": [],
            "message": ("DeepPrime ranking is only available in the Flask/server "
                        "build (needs the genet package). The full standalone "
                        "candidate list is shown below."),
        }

    return resp


# ── verify ──

def _build_variant(d):
    phase_map = {"0": 0, "1": 1, "2": 2}
    frame = d.get("frame", "start")
    var_phase = phase_map.get(frame)
    from_start = frame == "start"
    cn = str(d.get("codon_number", "")).strip()
    var_num = int(cn) if cn.isdigit() else None
    ex = _exon_kwargs(d)
    if d.get("mode", "pair") == "pair":
        return build_variant_from_pair(
            d["wt"], d["ed"], var_codon_phase=var_phase,
            wt_in_frame_from_start=from_start, var_codon_number=var_num, **ex)
    wt = d["wt"].strip().upper().replace(" ", "").replace("\n", "")
    pos, ref, alt = int(d["pos"]), d.get("ref", "").upper(), d.get("alt", "").upper()
    return build_variant_from_spec(
        wt, pos, ref, alt, var_codon_phase=var_phase,
        wt_in_frame_from_start=from_start, var_codon_number=var_num, **ex)


def run_verify(d: dict) -> dict:
    """Same contract as server.py /api/verify, minus Flask."""
    try:
        variant = _build_variant(d)
    except Exception as e:
        return {"ok": False, "error": "Could not build target: %s" % e}
    try:
        from silent_mutation.core.verify import verify_pegrna
        return verify_pegrna(
            variant,
            spacer=d.get("spacer", ""),
            extension=d.get("extension", ""),
            pbs_len=d.get("pbs_len_v", d.get("pbs_len")),
        )
    except Exception as e:
        traceback.print_exc()
        return {"ok": False, "error": "Verify failed: %s" % e}
