"""
silent_mutation.core.verify — verify a hand-designed pegRNA.

Given a Variant (target + 121 bp window) and a user's pegRNA written as
    spacer  +  3' extension (RTT + PBS, 5'->3')  +  PBS length
this decodes it to PAM / nick / PBS / RTT, extracts the silent substitutions
the user placed in the RTT, and validates each against the tool's own
synonymous-candidate set (find_silent_candidates) — the same ground truth the
main pipeline uses. No new biology: it reuses the pam_finder / types
conventions and the existing silent finder.

Why PBS length is required: from a concatenated extension alone the RTT|PBS
boundary is ambiguous (every WT-matching suffix length is consistent, and a
coincidental WT base at the RTT 3' end shifts it by +/-1). The designer always
knows their PBS length, so we take it explicitly and the split is exact.

Orientation (see pam_finder / types):
  seq_wt = CDS strand 5'->3', VAR_IDX = 60.
  spacer = protospacer, PAM strand 5'->3'.
    NGG: spacer == seq_wt[pam_pos-S : pam_pos], 'NGG' at [pam_pos:pam_pos+3]
    CCN: seq_wt[pam_pos+3 : pam_pos+3+S] == revcomp(spacer), 'CCN' at [pam_pos:pam_pos+3]
  3' extension as written 5'->3' (what DeepPrime / SynDesign output) is the
  REVERSE COMPLEMENT of the PAM-strand (PBS+RTT). i.e. rc(extension) == pbs_seq +
  rtt_seq (PAM strand), PBS first. So we reverse-complement the user's extension,
  split off pbs_len from the front as PBS (must equal WT pbs_seq), and the rest is
  the edited RTT. For robustness we also accept a user who pasted the forward
  PBS+RTT directly, and report which orientation matched ("orient").

Scope: substitution edits (ref_len == alt_len). Used by POST /api/verify.
Existing endpoints/behaviour are untouched.
"""

from __future__ import annotations

from typing import Optional

from silent_mutation.core.types import VAR_IDX, PamCandidate, Variant
from silent_mutation.core.silent_finder import find_silent_candidates
from silent_mutation.core.codon_utils import load_codon_table
from silent_mutation.io.genome_loader import reverse_complement

_COMP = str.maketrans("ACGT", "TGCA")


def _clean(s: str) -> str:
    return (s or "").upper().replace(" ", "").replace("\n", "").replace("\t", "")


# -- pure decode helpers (unit-testable without a Variant) --

def locate_geometry(seq_wt: str, spacer: str) -> list:
    """All (pam_pattern, pam_pos, nick_pos) in seq_wt consistent with `spacer`."""
    S = len(spacer)
    out = []
    start = 0
    while True:                                            # NGG (PAM strand == CDS)
        q = seq_wt.find(spacer, start)
        if q < 0:
            break
        pam_pos = q + S
        if pam_pos + 3 <= len(seq_wt) and seq_wt[pam_pos + 1:pam_pos + 3] == "GG":
            out.append(("NGG", pam_pos, pam_pos - 3))
        start = q + 1
    rc = reverse_complement(spacer)
    start = 0
    while True:                                            # CCN (PAM strand != CDS)
        r = seq_wt.find(rc, start)
        if r < 0:
            break
        pam_pos = r - 3
        if pam_pos >= 0 and seq_wt[pam_pos:pam_pos + 2] == "CC":
            out.append(("CCN", pam_pos, pam_pos + 6))
        start = r + 1
    return out


def wt_pbs_rtt(seq_wt: str, pattern: str, nick: int, pbs_len: int, rtt_len: int):
    """WT pbs_seq / rtt_seq (PAM-strand 5'->3') + rtt span, per pam_finder."""
    if pattern == "NGG":
        return seq_wt[nick - pbs_len:nick], seq_wt[nick:nick + rtt_len], nick, nick + rtt_len
    cds_pbs = seq_wt[nick:nick + pbs_len]
    cds_rtt = seq_wt[nick - rtt_len:nick]
    return reverse_complement(cds_pbs), reverse_complement(cds_rtt), nick - rtt_len, nick


def rtt_index_to_seqwt(pattern: str, j: int, rtt_start: int, rtt_end: int) -> int:
    return rtt_start + j if pattern == "NGG" else rtt_end - 1 - j


def seqwt_alt_from_rtt_base(pattern: str, base: str) -> str:
    return base if pattern == "NGG" else base.translate(_COMP)


# -- main entry --

def verify_pegrna(variant: Variant, spacer: str, extension: str, pbs_len: int,
                  codon_table: Optional[dict] = None) -> dict:
    """Decode + validate a hand-designed pegRNA. Returns a JSON-able verdict."""
    codon_table = codon_table or load_codon_table()
    seq_wt, seq_ed = variant.seq_wt, variant.seq_ed
    spacer = _clean(spacer)
    extension = _clean(extension)

    if not (16 <= len(spacer) <= 25):
        return {"ok": False, "error": "Spacer should be ~20 nt (the protospacer, 5'->3')."}
    if variant.ref_len != variant.alt_len:
        return {"ok": False, "error": "Verify currently supports substitution edits only."}
    try:
        pbs_len = int(pbs_len)
    except (TypeError, ValueError):
        return {"ok": False, "error": "Enter the PBS length (number of PBS nt in the extension)."}
    rtt_len = len(extension) - pbs_len
    if pbs_len < 1 or rtt_len < 1:
        return {"ok": False, "error": "PBS length must be between 1 and (extension length - 1)."}

    geoms = locate_geometry(seq_wt, spacer)
    if not geoms:
        return {"ok": False, "error": "Spacer not found next to an NGG/CCN PAM in the window. "
                "Check the spacer and that it is written 5'->3'."}

    ext_rc = reverse_complement(extension)   # rc(pegRNA 3' ext) == PBS+RTT (PAM strand)
    chosen = None
    for pattern, pam_pos, nick in geoms:
        if pattern == "NGG":
            if nick - pbs_len < 0 or nick + rtt_len > len(seq_wt):
                continue
        else:
            if nick - rtt_len < 0 or nick + pbs_len > len(seq_wt):
                continue
        pbs_wt, rtt_wt, rtt_s, rtt_e = wt_pbs_rtt(seq_wt, pattern, nick, pbs_len, rtt_len)
        if not (rtt_s <= VAR_IDX < rtt_e):
            continue
        je = (VAR_IDX - rtt_s) if pattern == "NGG" else (rtt_e - 1 - VAR_IDX)
        # forward (PAM strand) = PBS + RTT. Primary: user pasted the pegRNA 3' ext
        # (rev-comp, like DeepPrime) -> rc(extension) is forward. Fallback: user
        # pasted the forward PBS+RTT directly.
        for fwd, orient in ((ext_rc, "revcomp"), (extension, "forward")):
            pbs_user, rtt_user = fwd[:pbs_len], fwd[pbs_len:]
            pbs_match = (pbs_user == pbs_wt)
            edit_ok = (0 <= je < rtt_len) and \
                (seqwt_alt_from_rtt_base(pattern, rtt_user[je]) == variant.alt_allele)
            cand = dict(pattern=pattern, pam_pos=pam_pos, nick=nick, pbs_wt=pbs_wt, rtt_wt=rtt_wt,
                        rtt_s=rtt_s, rtt_e=rtt_e, rtt_user=rtt_user, pbs_match=pbs_match,
                        edit_ok=edit_ok, orient=orient)
            better = (chosen is None
                      or (pbs_match and edit_ok) and not (chosen["pbs_match"] and chosen["edit_ok"])
                      or (pbs_match and not chosen["pbs_match"]))
            if better:
                chosen = cand
            if pbs_match and edit_ok:
                break
        if chosen and chosen["pbs_match"] and chosen["edit_ok"]:
            break

    if chosen is None:
        return {"ok": False, "error": "The RTT at this PAM does not reach the edit with this "
                "PBS length - check the spacer / PBS length / RTT length."}

    pattern = chosen["pattern"]; nick = chosen["nick"]; pam_pos = chosen["pam_pos"]
    rtt_s, rtt_e = chosen["rtt_s"], chosen["rtt_e"]
    rtt_user, rtt_wt = chosen["rtt_user"], chosen["rtt_wt"]
    var_lo, var_hi = VAR_IDX, VAR_IDX + max(variant.ref_len, 1)

    user_silents = []
    for j in range(rtt_len):
        if rtt_user[j] == rtt_wt[j]:
            continue
        pos = rtt_index_to_seqwt(pattern, j, rtt_s, rtt_e)
        if var_lo <= pos < var_hi:
            continue
        user_silents.append((pos, seqwt_alt_from_rtt_base(pattern, rtt_user[j])))

    pam = PamCandidate(pam_pos=pam_pos, pam_pattern=pattern, nick_pos=nick,
                       pbs_seq=chosen["pbs_wt"], rtt_seq=rtt_wt,
                       rtt_start_in_seq=rtt_s, rtt_end_in_seq=rtt_e)
    valid = find_silent_candidates(variant=variant, pam=pam, codon_table=codon_table)
    valid_pos = {}
    for c in valid:
        for p, a in zip(c.silent_positions, c.silent_alt):
            valid_pos[p] = {"alt": a, "locale": c.locale, "priority": c.priority,
                            "codon": "%s\u2192%s" % (c.original_codon, c.mutated_codon), "aa": c.aa}

    checked, all_valid = [], True
    for pos, alt in user_silents:
        hit = valid_pos.get(pos)
        ok = bool(hit and hit["alt"] == alt)
        all_valid = all_valid and ok
        checked.append({
            "pos": pos, "alt": alt, "valid": ok,
            "locale": hit["locale"] if hit else "",
            "priority": hit["priority"] if hit else None,
            "codon": hit["codon"] if hit else "",
            "aa": hit["aa"] if hit else "",
            "reason": "" if ok else (
                "wrong base here (a synonymous change at this position would be \u2192%s)" % hit["alt"]
                if hit else
                "not a synonymous option here (likely non-synonymous, out of frame, or intron)"),
        })

    overall = chosen["edit_ok"] and chosen["pbs_match"] and all_valid and len(user_silents) > 0
    roles = _roles(seq_wt, pam, variant, [p for p, _ in user_silents])
    sil_off = [(p - rtt_s) if pattern == "NGG" else (rtt_e - 1 - p) for p, _ in user_silents]

    return {
        "ok": True, "overall": overall,
        "pattern": pattern, "strand": "+" if pattern == "NGG" else "-",
        "nick": nick, "pbs_len": pbs_len, "rtt_len": rtt_len,
        "pbs_seq": chosen["pbs_wt"], "rtt_seq_wt": rtt_wt, "rtt_user": rtt_user,
        "pbs_match": chosen["pbs_match"], "edit_ok": chosen["edit_ok"], "orient": chosen["orient"],
        "n_user_silents": len(user_silents), "silents": checked,
        "n_valid_options": len(valid_pos),
        "seq_wt": seq_wt, "roles": roles, "sil_off": sil_off, "rtt": rtt_user, "nick_pos": nick,
    }


def _roles(seq_wt: str, pam, variant: Variant, silent_positions: list) -> list:
    L = len(seq_wt)
    track = ["bg"] * L
    rtt_s, rtt_e = pam.rtt_start_in_seq, pam.rtt_end_in_seq
    pam_s, pam_e = pam.pam_pos, pam.pam_pos + 3
    if pam.pam_pattern == "NGG":
        pbs_s, pbs_e = pam.nick_pos - len(pam.pbs_seq), pam.nick_pos
    else:
        pbs_s, pbs_e = pam.nick_pos, pam.nick_pos + len(pam.pbs_seq)
    for i in range(max(rtt_s, 0), min(rtt_e, L)):
        track[i] = "rtt"
    for i in range(max(pbs_s, 0), min(pbs_e, L)):
        track[i] = "pbs"
    for i in range(max(pam_s, 0), min(pam_e, L)):
        track[i] = "pam"
    for i in range(VAR_IDX, min(VAR_IDX + max(variant.ref_len, 1), L)):
        track[i] = "variant"
    for p in silent_positions:
        if 0 <= p < L:
            track[p] = "silent"
    return track
