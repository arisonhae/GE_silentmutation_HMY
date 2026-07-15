"""
Silent (synonymous) mutation candidate finder.

For a given (Variant, PamCandidate) pair, scan the RTT region for codons
where one or more nucleotide substitutions would yield a synonymous codon
(same amino acid). These candidates serve as sequence markers and as PAM
disruptors.

Codons are read from Variant.codon_lookup (position-keyed), not by slicing
seq_wt; positions absent from the lookup are skipped. See sequence_loader for
how the lookup is built and what it omits.

Algorithm
---------
1. For each position in RTT ∩ seq_wt:
   - Convert seq_wt index → genomic position
   - Look up the codon via variant.codon_lookup
   - If None (intron base, UTR), skip
2. Deduplicate by CodonInfo.codon_index
3. Skip codons overlapping the variant
4. For each synonymous alternative:
   - Compute silent base positions in genome → seq_wt
   - Verify all silent bases are in seq_wt AND in the RTT region
   - Verify none of them are intron bases (codon_lookup contains them)
   - Classify locale (PAM / LHA / RHA) and emit a SilentCandidate

Coordinates everywhere are 0-based half-open into Variant.seq_wt (CDS
strand 5'->3'), matching types.py and pam_finder.py.
"""

from __future__ import annotations

from typing import Optional

from silent_mutation.core.codon_utils import load_codon_table
from silent_mutation.core.types import (
    PamCandidate,
    SilentCandidate,
    Variant,
)


# ─────────────────────────────────────────────────────────────────────────────
# Helpers
# ─────────────────────────────────────────────────────────────────────────────

def _variant_span(var_idx: int, var_ref_len: int) -> tuple[int, int]:
    """Inclusive-exclusive variant span on seq_wt. ref_len=0 → single pos."""
    if var_ref_len == 0:
        return var_idx, var_idx + 1
    return var_idx, var_idx + var_ref_len


def _codon_overlaps_variant_by_genome(
    codon_genomic_positions: tuple[int, int, int],
    variant_genome_pos: int,
    variant_ref_len: int,
) -> bool:
    """
    Does the codon (given by its 3 genomic positions) overlap the variant's
    genomic span? Works regardless of strand since we compare on the plus-
    strand genomic coordinate axis.

    For insertions (ref_len=0), variant occupies a single genomic position.
    """
    if variant_ref_len == 0:
        var_lo, var_hi = variant_genome_pos, variant_genome_pos + 1
    else:
        var_lo, var_hi = variant_genome_pos, variant_genome_pos + variant_ref_len
    return any(var_lo <= g < var_hi for g in codon_genomic_positions)


def _classify_locale(
    pos_in_seq: int,
    pam: PamCandidate,
    var_idx: int,
    var_ref_len: int,
) -> str:
    """
    Classify a seq_wt position as 'PAM' / 'LHA' / 'RHA' relative to the
    PAM and variant. Symmetric for NGG/CCN — uses pam_pos and variant
    location, not raw direction.
    """
    span_lo, span_hi = _variant_span(var_idx, var_ref_len)

    # PAM triplet: [pam_pos, pam_pos + 3)
    if pam.pam_pos <= pos_in_seq < pam.pam_pos + 3:
        return "PAM"

    if pam.pam_pattern == "NGG":
        # PAM strand 5'->3' == CDS strand 5'->3'
        # Along PAM strand: nick → variant → RTT end
        if pos_in_seq < span_lo:
            return "LHA"
        if pos_in_seq >= span_hi:
            return "RHA"
        return "RHA"

    # CCN: PAM strand antiparallel to CDS strand.
    # Along PAM strand 5'->3', higher seq_wt index comes first (closer to nick),
    # then variant, then lower indices.
    if pos_in_seq >= span_hi:
        return "LHA"
    if pos_in_seq < span_lo:
        return "RHA"
    return "RHA"


_LOCALE_PRIORITY = {"PAM": 1, "LHA": 2, "RHA": 3}


def _enumerate_synonymous_codons(
    original: str,
    table: dict[str, str],
) -> list[str]:
    """All codons synonymous to `original`, EXCLUDING `original` itself."""
    target_aa = table.get(original)
    if target_aa is None:
        return []
    return [c for c, aa in table.items() if aa == target_aa and c != original]


def _diff_positions(a: str, b: str) -> list[int]:
    """Indices (0..2) where two 3-letter codons differ."""
    return [i for i in range(3) if a[i] != b[i]]


# ─────────────────────────────────────────────────────────────────────────────
# Main API
# ─────────────────────────────────────────────────────────────────────────────

def find_silent_candidates(
    variant: Variant,
    pam: PamCandidate,
    codon_table: Optional[dict[str, str]] = None,
) -> list[SilentCandidate]:
    """
    Find all silent-mutation candidates inside the RTT region of one PAM.

    Returns
    -------
    list[SilentCandidate]
        One entry per (codon, synonymous-alternative) pair satisfying:
          - codon does NOT overlap the variant (genomic-coordinate check)
          - every silent edit position is a CDS base (in codon_lookup)
          - every silent edit position lies inside the RTT region
        Sorted by (priority asc, leftmost silent_position asc).
    """
    if codon_table is None:
        codon_table = load_codon_table()

    if not variant.codon_lookup:
        return []

    candidates: list[SilentCandidate] = []
    seen_codon_indices: set[int] = set()

    rtt_start = pam.rtt_start_in_seq
    rtt_end = pam.rtt_end_in_seq

    # Walk every position in the RTT region. For each, find which codon
    # (if any) contains the base at that seq_wt position. Dedupe by
    # codon_index so each codon is examined once.
    for seq_idx in range(rtt_start, rtt_end):
        gpos = variant.seq_idx_to_genome_pos(seq_idx)
        codon_info = variant.codon_lookup.get(gpos)
        if codon_info is None:
            # Intron base, UTR, or outside CDS — not a codon base
            continue
        if codon_info.codon_index in seen_codon_indices:
            continue
        seen_codon_indices.add(codon_info.codon_index)

        # Skip codons overlapping the variant. Use genomic-coordinate check
        # because splice-junction codons have non-contiguous seq_wt indices.
        if _codon_overlaps_variant_by_genome(
            codon_info.codon_genomic_positions,
            variant.genome_pos,
            variant.ref_len,
        ):
            continue

        # Enumerate synonymous alternatives
        for alt_codon in _enumerate_synonymous_codons(codon_info.codon, codon_table):
            diff_idx = _diff_positions(codon_info.codon, alt_codon)

            # For each differing position, find the seq_wt index of that base.
            # codon_info.codon_genomic_positions gives the 3 genomic positions
            # of the codon in mRNA 5'->3' order (frame 0, 1, 2).
            silent_seq_indices: list[int] = []
            in_rtt = True
            for fi in diff_idx:
                gp = codon_info.codon_genomic_positions[fi]
                si = variant.genome_pos_to_seq_idx(gp)
                if si is None:
                    # silent base falls outside seq_wt window — can't edit it
                    in_rtt = False
                    break
                if not (rtt_start <= si < rtt_end):
                    # silent base is in seq_wt but outside RTT — can't edit it
                    in_rtt = False
                    break
                silent_seq_indices.append(si)
            if not in_rtt:
                continue

            silent_ref = "".join(codon_info.codon[i] for i in diff_idx)
            silent_alt = "".join(alt_codon[i] for i in diff_idx)

            # Locale: take the strongest classification among the silent bases.
            # PAM > LHA > RHA.
            locales = [
                _classify_locale(si, pam, variant.var_idx, variant.ref_len)
                for si in silent_seq_indices
            ]
            if "PAM" in locales:
                locale = "PAM"
            elif "LHA" in locales:
                locale = "LHA"
            else:
                locale = "RHA"
            priority = _LOCALE_PRIORITY[locale]

            candidates.append(SilentCandidate(
                pam=pam,
                silent_positions=silent_seq_indices,
                silent_ref=silent_ref,
                silent_alt=silent_alt,
                original_codon=codon_info.codon,
                mutated_codon=alt_codon,
                aa=codon_info.aa,
                locale=locale,           # type: ignore[arg-type]
                priority=priority,
            ))

    # Stable sort: priority asc, then leftmost silent position asc
    candidates.sort(key=lambda c: (c.priority, min(c.silent_positions)))
    return candidates
