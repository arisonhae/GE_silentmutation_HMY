"""
silent_mutation.webtool.server — Flask backend.

Serves a hand-built HTML front-end (index.html) and JSON endpoints that run
the silent-mutation pipeline.

Run (in a 'genet' conda env with flask installed), from the repo root:
    PYTHONPATH=. python silent_mutation/webtool/server.py
Then expose, e.g.:
    ngrok http 8502

Endpoints:
    GET  /              -> index.html (the web tool)
    POST /api/analyze   -> {variant, candidates[], deepprime} as JSON
    POST /api/verify    -> verdict for a hand-designed pegRNA, as JSON

Exon-aware mode
---------------
When the front-end's reading-frame bar is in exon-aware mode it sends
exon_start / exon_end / codon_start (0-based indices into the WT string) and an
optional cds_coord. These flow into build_variant_from_pair / build_variant_from_spec
AND into run_deepprime_silent, so the codon lookup is built over the exon only
everywhere — an intronic edit carries its silent markers in the adjacent exon,
in both the standalone candidate list and the DeepPrime ranking.
"""

from __future__ import annotations

import os
import traceback

from flask import Flask, request, jsonify, send_from_directory

from silent_mutation.io.sequence_loader import (
    build_variant_from_pair, build_variant_from_spec,
)
from silent_mutation.core.pegrna_builder import build_outputs_for_variant
from silent_mutation.io.genome_loader import reverse_complement
from silent_mutation.io.deepprime_runner import (
    run_deepprime_silent, DeepPrimeUnavailable,
)

HERE = os.path.dirname(os.path.abspath(__file__))
app = Flask(__name__)


# ── visualization: per-base roles (silent > variant > PAM > PBS > RTT > bg) ──

def role_track(row) -> list[str]:
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
    """Return (pbsrtt_forward, ext_rc).
    pbsrtt_forward = pbs_seq + edited_rtt   (PAM strand 5'->3', = DeepPrime 'pbsrtt')
    ext_rc         = reverse_complement(pbsrtt_forward)  (= DeepPrime 'Extension Top',
                     the pegRNA 3' extension as written 5'->3'). Edited RTT carries
                     variant + silent."""
    try:
        if len(row.seq_ed_with_silent) != len(row.seq_wt):
            return "", ""                              # indel: coords shift; skip for now
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
    """Pull exon-aware args from the request body, if the bar sent them.

    Returns {} for the normal contiguous path. All three of exon_start /
    exon_end / codon_start must be present together; cds_coord is optional.
    """
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


# ── routes ──

@app.route("/")
def index():
    return send_from_directory(HERE, "index.html")


@app.route("/api/analyze", methods=["POST"])
def analyze():
    d = request.get_json(force=True)
    try:
        phase_map = {"0": 0, "1": 1, "2": 2}
        frame = d.get("frame", "start")
        var_phase = phase_map.get(frame)            # None => infer from boundary
        from_start = frame == "start"
        cn = str(d.get("codon_number", "")).strip()
        var_num = int(cn) if cn.isdigit() else None

        ex = _exon_kwargs(d)                         # {} unless bar is in exon mode

        mode = d.get("mode", "pair")
        if mode == "pair":
            wt, ed = d["wt"], d["ed"]
            variant = build_variant_from_pair(
                wt, ed, var_codon_phase=var_phase,
                wt_in_frame_from_start=from_start, var_codon_number=var_num,
                **ex,
            )
            dp_wt, dp_ed = wt, ed
        else:
            wt = d["wt"].strip().upper().replace(" ", "").replace("\n", "")
            pos, ref, alt = int(d["pos"]), d.get("ref", "").upper(), d.get("alt", "").upper()
            variant = build_variant_from_spec(
                wt, pos, ref, alt, var_codon_phase=var_phase,
                wt_in_frame_from_start=from_start, var_codon_number=var_num,
                **ex,
            )
            dp_wt = wt
            dp_ed = wt[:pos] + alt + wt[pos + len(ref):]
    except Exception as e:
        return jsonify({"ok": False, "error": str(e)}), 400

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

    # DeepPrime ranking (genet). Works in both contiguous and exon-aware mode:
    # the exon args flow into the variant build, so silent markers per pegRNA come
    # from the exon lookup just like the standalone path. DeepPrime scoring itself
    # is frame-independent (it works on the sequence).
    if d.get("use_dp"):
        dp = d.get("dp", {})
        eff_phase = info.frame if info else 0
        try:
            pegrnas = run_deepprime_silent(
                dp_wt, dp_ed,
                var_codon_phase=eff_phase, wt_in_frame_from_start=False,
                pe_system=dp.get("pe_system", "PE2max"),
                cell_type=dp.get("cell_type", "HEK293T"),
                rtt_max=max_rtt, top_n=int(dp.get("top_n", 20)),
                dp_path=dp.get("path", ""),
                **ex,
            )
            ser = []
            for p in pegrnas:
                sp = {k: v for k, v in p.items() if k != "outputs"}
                sp["silents"] = [candidate_dict(o, i + 1) for i, o in enumerate(p["outputs"])]
                ser.append(sp)
            resp["deepprime"] = {
                "available": True,
                "pegrnas": ser,
                "message": "" if ser else "No pegRNAs found (no NGG PAM in range).",
            }
        except DeepPrimeUnavailable as e:
            resp["deepprime"] = {"available": False, "pegrnas": [], "message": str(e)}
        except Exception as e:
            resp["deepprime"] = {"available": False, "pegrnas": [],
                                 "message": "DeepPrime run failed: %s" % e}

    return jsonify(resp)


def _build_variant(d):
    """Build a Variant from the request body — same contract as /api/analyze.

    Used by /api/verify so the verification target is constructed identically
    (pair/spec mode, frame/codon_number, exon-aware args). Does not alter the
    existing /api/analyze path.
    """
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


@app.route("/api/verify", methods=["POST"])
def verify():
    """Verify a hand-designed pegRNA (spacer + 3' extension + PBS length)
    against the same WT/edit window. Independent of /api/analyze."""
    d = request.get_json(force=True)
    try:
        variant = _build_variant(d)
    except Exception as e:
        return jsonify({"ok": False, "error": "Could not build target: %s" % e}), 400
    try:
        from silent_mutation.core.verify import verify_pegrna
        res = verify_pegrna(
            variant,
            spacer=d.get("spacer", ""),
            extension=d.get("extension", ""),
            pbs_len=d.get("pbs_len_v", d.get("pbs_len")),
        )
        return jsonify(res)
    except Exception as e:
        traceback.print_exc()
        return jsonify({"ok": False, "error": "Verify failed: %s" % e})


if __name__ == "__main__":
    app.run(host="0.0.0.0", port=8502, debug=False)
