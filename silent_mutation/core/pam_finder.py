"""
PAM finder for prime editing.

Scans seq_wt for NGG (CDS strand) and CCN (PAM strand) PAM sites near the
variant, computes nick / PBS / RTT coordinates per types.py conventions,
and returns PamCandidate objects.

Coordinate conventions (recap from types.py)
--------------------------------------------
All positions are 0-based half-open indices into Variant.seq_wt.
seq_wt is stored 5'->3' on the CDS strand.

For pam_pattern == 'NGG':
    PAM strand == CDS strand
    pam_pos    = index of the N in NGG
    nick_pos   = pam_pos - 3
    RTT region = [nick_pos, nick_pos + rtt_len)         (PAM included)
    PBS region = [nick_pos - pbs_len, nick_pos)
    pbs_seq, rtt_seq taken directly from seq_wt

For pam_pattern == 'CCN':
    PAM strand != CDS strand
    pam_pos    = index of the first C in CCN
    nick_pos   = pam_pos + 6
    RTT region = [nick_pos - rtt_len, nick_pos)         (PAM included)
    PBS region = [nick_pos, nick_pos + pbs_len)
    pbs_seq, rtt_seq are the reverse complement of the seq_wt slice
        (so they are stored 5'->3' on the PAM strand, matching DeepPrime
        convention)

RTT length policy
-----------------
SynDesign-style dynamic length: each PAM gets the largest RTT that fits
within seq_wt (capped by max_rtt_len, default 40). This means any silent
candidate that would be reachable at any RTT length <= max_rtt_len is
captured by a single scan at max — shorter RTTs are strict subsets of
the same window starting at nick_pos.

Variant-in-RTT filter
---------------------
By default, PAM candidates whose RTT region does not cover the variant
are dropped (they cannot install the intended edit). The variant occupies
[var_idx, var_idx + var_ref_len) in seq_wt; for an insertion (ref_len=0)
we treat the insertion site as the single position var_idx.
"""

from __future__ import annotations

import re
from typing import Optional

from silent_mutation.core.types import (
    FLANK,
    VAR_IDX,
    PamCandidate,
    Variant,
)
from silent_mutation.io.genome_loader import reverse_complement

DEFAULT_PBS_LEN = 13
DEFAULT_MAX_RTT_LEN = 40

# PAM patterns. We use simple regex with explicit base classes so that
# 'N' positions only match real ACGT (not literal 'N' / masked bases).
_NGG_RE = re.compile(r"(?=([ACGT]GG))")
_CCN_RE = re.compile(r"(?=(CC[ACGT]))")


# ─────────────────────────────────────────────────────────────────────────────
# Helpers
# ─────────────────────────────────────────────────────────────────────────────

def _variant_span(var_idx: int, var_ref_len: int) -> tuple[int, int]:
    """
    Return [span_lo, span_hi) — the inclusive range of seq_wt positions
    that the variant occupies. For an insertion (ref_len=0) this collapses
    to a single position [var_idx, var_idx+1) so the "is the variant in
    the RTT" check still does something sensible.
    """
    if var_ref_len == 0:
        return var_idx, var_idx + 1
    return var_idx, var_idx + var_ref_len


def _variant_in_range(
    range_start: int, range_end: int,
    var_idx: int, var_ref_len: int,
) -> bool:
    """True iff the variant span overlaps [range_start, range_end)."""
    span_lo, span_hi = _variant_span(var_idx, var_ref_len)
    return span_lo < range_end and span_hi > range_start


# ─────────────────────────────────────────────────────────────────────────────
# Per-strand scan
# ─────────────────────────────────────────────────────────────────────────────

def _scan_ngg(
    seq: str,
    var_idx: int,
    var_ref_len: int,
    pbs_len: int,
    max_rtt_len: int,
    require_variant_in_rtt: bool,
) -> list[PamCandidate]:
    """Find all NGG PAM candidates on the CDS strand."""
    out: list[PamCandidate] = []
    seq_len = len(seq)

    for m in _NGG_RE.finditer(seq):
        pam_pos = m.start()

        # nick is 3bp upstream of N (5' side on CDS strand)
        nick_pos = pam_pos - 3

        # PBS sits to the LEFT of nick on CDS strand
        pbs_start = nick_pos - pbs_len
        pbs_end = nick_pos

        # RTT sits to the RIGHT of nick (PAM included)
        rtt_start = nick_pos
        # Dynamic length: bounded by seq_wt edge
        rtt_len = min(max_rtt_len, seq_len - rtt_start)
        rtt_end = rtt_start + rtt_len

        # Bounds checks
        if pbs_start < 0:
            continue
        if rtt_len <= 0:
            continue
        # RTT must at minimum extend past the PAM (3bp); otherwise the
        # candidate is degenerate.
        if rtt_len < 3:
            continue

        # Variant-in-RTT filter
        if require_variant_in_rtt and not _variant_in_range(
            rtt_start, rtt_end, var_idx, var_ref_len
        ):
            continue

        pbs_seq = seq[pbs_start:pbs_end]
        rtt_seq = seq[rtt_start:rtt_end]

        # Defensive: skip if any N or non-ACGT slipped in
        if not _is_clean_acgt(pbs_seq) or not _is_clean_acgt(rtt_seq):
            continue

        out.append(PamCandidate(
            pam_pos=pam_pos,
            pam_pattern="NGG",
            nick_pos=nick_pos,
            pbs_seq=pbs_seq,
            rtt_seq=rtt_seq,
            rtt_start_in_seq=rtt_start,
            rtt_end_in_seq=rtt_end,
        ))

    return out


def _scan_ccn(
    seq: str,
    var_idx: int,
    var_ref_len: int,
    pbs_len: int,
    max_rtt_len: int,
    require_variant_in_rtt: bool,
) -> list[PamCandidate]:
    """Find all CCN PAM candidates (i.e. NGG on the opposite strand)."""
    out: list[PamCandidate] = []
    seq_len = len(seq)

    for m in _CCN_RE.finditer(seq):
        pam_pos = m.start()       # first C of CCN

        # nick is 6bp downstream of pam_pos on CDS strand index axis
        # (corresponds to "3bp 5' of PAM on PAM strand", same biology
        # as the NGG case — see types.py docstring).
        nick_pos = pam_pos + 6

        # On CDS strand: RTT is to the LEFT of nick (and includes the CCN);
        # PBS is to the RIGHT of nick.
        rtt_end = nick_pos
        rtt_len = min(max_rtt_len, rtt_end)        # bounded by left edge
        rtt_start = rtt_end - rtt_len

        pbs_start = nick_pos
        pbs_end = nick_pos + pbs_len

        # Bounds checks
        if pbs_end > seq_len:
            continue
        if rtt_len <= 0:
            continue
        if rtt_len < 3:
            continue

        # Variant-in-RTT filter
        if require_variant_in_rtt and not _variant_in_range(
            rtt_start, rtt_end, var_idx, var_ref_len
        ):
            continue

        # Slice on CDS strand, then RC to put pbs/rtt in PAM-strand 5'->3'
        cds_pbs = seq[pbs_start:pbs_end]
        cds_rtt = seq[rtt_start:rtt_end]
        if not _is_clean_acgt(cds_pbs) or not _is_clean_acgt(cds_rtt):
            continue

        pbs_seq = reverse_complement(cds_pbs)
        rtt_seq = reverse_complement(cds_rtt)

        out.append(PamCandidate(
            pam_pos=pam_pos,
            pam_pattern="CCN",
            nick_pos=nick_pos,
            pbs_seq=pbs_seq,
            rtt_seq=rtt_seq,
            rtt_start_in_seq=rtt_start,
            rtt_end_in_seq=rtt_end,
        ))

    return out


def _is_clean_acgt(s: str) -> bool:
    """True if `s` consists entirely of A/C/G/T."""
    return bool(s) and set(s).issubset("ACGT")


# ─────────────────────────────────────────────────────────────────────────────
# Public API
# ─────────────────────────────────────────────────────────────────────────────

def find_pam_candidates(
    seq: str,
    var_idx: int = VAR_IDX,
    var_ref_len: int = 1,
    pbs_len: int = DEFAULT_PBS_LEN,
    max_rtt_len: int = DEFAULT_MAX_RTT_LEN,
    require_variant_in_rtt: bool = True,
) -> list[PamCandidate]:
    """
    Scan a CDS-strand sequence for NGG and CCN PAM sites near the variant
    and return all PamCandidate objects whose RTT region (with the chosen
    max_rtt_len) covers the variant.

    Parameters
    ----------
    seq : str
        seq_wt in CDS-strand 5'->3' orientation. Must be uppercase ACGT.
    var_idx : int
        0-based index of the variant's first base in seq.
    var_ref_len : int
        REF allele length on the CDS strand. 0 for pure insertion.
    pbs_len : int
        PBS length to use for the candidate (does NOT affect silent finding;
        included so PamCandidate carries a usable pbs_seq for downstream
        DeepPrime input). Default 13 (DeepPrime mid-range default).
    max_rtt_len : int
        Upper bound on RTT length. Each PAM gets the largest RTT that fits
        within seq, up to this cap. Default 40 (SynDesign standard).
    require_variant_in_rtt : bool
        If True (default), PAMs whose RTT does not span the variant are
        excluded — they cannot install the edit. Setting False returns all
        PAMs unfiltered (useful for diagnostics).

    Returns
    -------
    list[PamCandidate], in order of increasing pam_pos. Both NGG and CCN
    candidates are interleaved by position.
    """
    if not seq:
        return []

    seq_upper = seq.upper()

    ngg = _scan_ngg(
        seq_upper, var_idx, var_ref_len, pbs_len, max_rtt_len,
        require_variant_in_rtt,
    )
    ccn = _scan_ccn(
        seq_upper, var_idx, var_ref_len, pbs_len, max_rtt_len,
        require_variant_in_rtt,
    )

    # Merge and sort by pam_pos for stable downstream iteration
    merged = ngg + ccn
    merged.sort(key=lambda c: (c.pam_pos, c.pam_pattern))
    return merged


def find_pam_candidates_for_variant(
    variant: Variant,
    pbs_len: int = DEFAULT_PBS_LEN,
    max_rtt_len: int = DEFAULT_MAX_RTT_LEN,
    require_variant_in_rtt: bool = True,
) -> list[PamCandidate]:
    """
    Convenience wrapper: extract var_idx and ref_len from a Variant.

    Note: scans seq_wt (the unedited window). The PAM/protospacer recognition
    happens before editing, so seq_wt is the biologically correct source.
    """
    return find_pam_candidates(
        seq=variant.seq_wt,
        var_idx=variant.var_idx,
        var_ref_len=variant.ref_len,
        pbs_len=pbs_len,
        max_rtt_len=max_rtt_len,
        require_variant_in_rtt=require_variant_in_rtt,
    )
