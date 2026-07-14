"""
pegRNA output builder (v2).

Combines (Variant, PamCandidate, SilentCandidate) triples into PegRNAOutput
rows. Updated for the v2 silent_finder API and to annotate each output row
with the variant-induced protein change (e.g. D23Y).
"""

from __future__ import annotations

from typing import Optional


from silent_mutation.core.codon_utils import load_codon_table
from silent_mutation.core.types import (
    VAR_IDX,
    PamCandidate,
    PegRNAOutput,
    SilentCandidate,
    Variant,
)
from silent_mutation.core.pam_finder import find_pam_candidates_for_variant
from silent_mutation.core.silent_finder import find_silent_candidates


# ─────────────────────────────────────────────────────────────────────────────
# Sequence application
# ─────────────────────────────────────────────────────────────────────────────

def _apply_substitutions(seq: str, positions: list[int], alts: str) -> str:
    if len(positions) != len(alts):
        raise ValueError(
            f"positions ({len(positions)}) and alts ({len(alts)}) length mismatch"
        )
    chars = list(seq)
    for p, a in zip(positions, alts):
        if not (0 <= p < len(chars)):
            raise IndexError(f"silent position {p} out of bounds (len={len(chars)})")
        chars[p] = a
    return "".join(chars)


def _seq_wt_pos_to_seq_ed_pos(pos: int, ref_len: int, alt_len: int) -> int:
    var_lo = VAR_IDX
    var_hi_wt = VAR_IDX + ref_len
    if var_lo <= pos < var_hi_wt:
        raise ValueError(
            f"seq_wt position {pos} lies inside variant span "
            f"[{var_lo}, {var_hi_wt}); cannot map to seq_ed."
        )
    if pos < var_lo:
        return pos
    return pos + (alt_len - ref_len)


def apply_silent_to_seq_ed(
    seq_ed: str,
    silent_positions: list[int],
    silent_alt: str,
    ref_len: int,
    alt_len: int,
) -> str:
    ed_positions = [
        _seq_wt_pos_to_seq_ed_pos(p, ref_len, alt_len) for p in silent_positions
    ]
    return _apply_substitutions(seq_ed, ed_positions, silent_alt)


# ─────────────────────────────────────────────────────────────────────────────
# Variant protein change
# ─────────────────────────────────────────────────────────────────────────────

def compute_variant_protein_change(
    variant: Variant,
    codon_table: Optional[dict[str, str]] = None,
) -> dict:
    """
    Compute the protein change caused by the variant at its own codon.

    Works for substitution variants (1-3bp) within a single codon. For
    indels or multi-codon-spanning substitutions, the simple "single
    codon → single AA" model doesn't apply cleanly; in that case we
    return blanks and let the caller decide how to label it.

    Returns a dict with keys:
        variant_protein_change : str  (e.g. "D23Y", or "" if N/A)
        variant_wt_codon       : str
        variant_mut_codon      : str
        variant_wt_aa          : str
        variant_mut_aa         : str
        variant_codon_index    : int  (0 if N/A)
    """
    blank = {
        "variant_protein_change": "",
        "variant_wt_codon": "",
        "variant_mut_codon": "",
        "variant_wt_aa": "",
        "variant_mut_aa": "",
        "variant_codon_index": 0,
    }

    if codon_table is None:
        codon_table = load_codon_table()

    info = variant.codon_at_variant()
    if info is None:
        return blank

    # Only handle simple substitutions confined to one codon
    if variant.ref_len != variant.alt_len or variant.ref_len == 0:
        # indel — skip protein change annotation
        return {**blank,
                "variant_wt_codon": info.codon,
                "variant_wt_aa": info.aa,
                "variant_codon_index": info.codon_index}

    # Substitution: check all variant bases fall in the same codon.
    # Each CDS-strand variant base (offset 0..ref_len-1) sits at seq_wt index
    # VAR_IDX+offset. seq_idx_to_genome_pos handles strand AND the leftmost-
    # anchor (ref_len-1) correction for minus-strand multi-bp variants, so we
    # no longer hand-roll genome_pos ± offset (which was wrong for minus-strand
    # 2-3bp substitutions).
    same_codon = True
    for offset in range(variant.ref_len):
        gpos = variant.seq_idx_to_genome_pos(VAR_IDX + offset)
        other_info = variant.codon_lookup.get(gpos)
        if other_info is None or other_info.codon_index != info.codon_index:
            same_codon = False
            break

    if not same_codon:
        # Variant spans codon boundary — would change 2 AAs
        return {**blank,
                "variant_wt_codon": info.codon,
                "variant_wt_aa": info.aa,
                "variant_codon_index": info.codon_index}

    # Build mut_codon by replacing the appropriate frame positions
    mut_codon_chars = list(info.codon)
    for offset in range(variant.ref_len):
        gpos = variant.seq_idx_to_genome_pos(VAR_IDX + offset)
        other_info = variant.codon_lookup[gpos]
        mut_codon_chars[other_info.frame] = variant.alt_allele[offset]
    mut_codon = "".join(mut_codon_chars)
    mut_aa = codon_table.get(mut_codon, "?")

    protein_change = f"{info.aa}{info.codon_index}{mut_aa}"

    return {
        "variant_protein_change": protein_change,
        "variant_wt_codon": info.codon,
        "variant_mut_codon": mut_codon,
        "variant_wt_aa": info.aa,
        "variant_mut_aa": mut_aa,
        "variant_codon_index": info.codon_index,
    }


# ─────────────────────────────────────────────────────────────────────────────
# Single-row builder
# ─────────────────────────────────────────────────────────────────────────────



def _compute_rtt_splice_context(variant: Variant, pam: PamCandidate) -> dict:
    """
    Walk every position in the RTT region and classify it using
    variant.codon_lookup. Returns counts useful for understanding why
    a given pegRNA has few/many silent candidates.

    Returns
    -------
    {
        "cds_bases_in_rtt": int,
        "rtt_covers_intron": bool,
        "splice_junction_codons_in_rtt": int,
    }
    """
    rtt_start = pam.rtt_start_in_seq
    rtt_end = pam.rtt_end_in_seq
    rtt_len = rtt_end - rtt_start

    if not variant.codon_lookup:
        return {
            "cds_bases_in_rtt": 0,
            "rtt_covers_intron": True,
            "splice_junction_codons_in_rtt": 0,
        }

    cds_count = 0
    junction_codon_indices = set()
    seen_codon_indices = set()

    for seq_idx in range(rtt_start, rtt_end):
        gpos = variant.seq_idx_to_genome_pos(seq_idx)
        info = variant.codon_lookup.get(gpos)
        if info is None:
            continue  # intron base, UTR — not CDS
        cds_count += 1
        if info.codon_index in seen_codon_indices:
            continue
        seen_codon_indices.add(info.codon_index)
        # Splice-junction codon? Check whether the three genomic positions
        # are contiguous.
        positions = info.codon_genomic_positions
        if variant.cds_strand == "+":
            contiguous = (positions[1] == positions[0] + 1 and
                          positions[2] == positions[1] + 1)
        else:
            contiguous = (positions[1] == positions[0] - 1 and
                          positions[2] == positions[1] - 1)
        if not contiguous:
            junction_codon_indices.add(info.codon_index)

    return {
        "cds_bases_in_rtt": cds_count,
        "rtt_covers_intron": (cds_count < rtt_len),
        "splice_junction_codons_in_rtt": len(junction_codon_indices),
    }


def build_pegrna_output(
    variant: Variant,
    pam: PamCandidate,
    silent: SilentCandidate,
    variant_protein: Optional[dict] = None,
    splice_ctx: Optional[dict] = None,
) -> PegRNAOutput:
    """
    Assemble a single PegRNAOutput row. `variant_protein` and `splice_ctx`
    may be passed pre-computed to avoid recomputing per (pam, silent) combo;
    fall back to computing here if None.
    """
    if splice_ctx is None:
        splice_ctx = _compute_rtt_splice_context(variant, pam)

    seq_ed_with_silent = apply_silent_to_seq_ed(
        seq_ed=variant.seq_ed,
        silent_positions=silent.silent_positions,
        silent_alt=silent.silent_alt,
        ref_len=variant.ref_len,
        alt_len=variant.alt_len,
    )

    if variant_protein is None:
        variant_protein = compute_variant_protein_change(variant)

    return PegRNAOutput(
        variant_id=variant.variant_id,
        gene_symbol=variant.gene_symbol,
        transcript_id=variant.transcript_id,
        chrom=variant.chrom,
        genome_pos=variant.genome_pos,
        ref_allele=variant.ref_allele,
        alt_allele=variant.alt_allele,
        cds_strand=variant.cds_strand,
        variant_type=variant.variant_type,
        ref_len=variant.ref_len,
        alt_len=variant.alt_len,
        seq_wt=variant.seq_wt,
        seq_ed=variant.seq_ed,
        seq_ed_with_silent=seq_ed_with_silent,
        pam_pos=pam.pam_pos,
        pam_pattern=pam.pam_pattern,
        nick_pos=pam.nick_pos,
        pbs_seq=pam.pbs_seq,
        rtt_seq_wt=pam.rtt_seq,
        rtt_start_in_seq=pam.rtt_start_in_seq,
        rtt_end_in_seq=pam.rtt_end_in_seq,
        silent_positions=list(silent.silent_positions),
        silent_ref=silent.silent_ref,
        silent_alt=silent.silent_alt,
        original_codon=silent.original_codon,
        mutated_codon=silent.mutated_codon,
        aa=silent.aa,
        locale=silent.locale,
        priority=silent.priority,
        variant_protein_change=variant_protein["variant_protein_change"],
        variant_wt_codon=variant_protein["variant_wt_codon"],
        variant_mut_codon=variant_protein["variant_mut_codon"],
        variant_wt_aa=variant_protein["variant_wt_aa"],
        variant_mut_aa=variant_protein["variant_mut_aa"],
        variant_codon_index=variant_protein["variant_codon_index"],
        cds_bases_in_rtt=splice_ctx["cds_bases_in_rtt"],
        rtt_covers_intron=splice_ctx["rtt_covers_intron"],
        splice_junction_codons_in_rtt=splice_ctx["splice_junction_codons_in_rtt"],
    )


# ─────────────────────────────────────────────────────────────────────────────
# Per-variant orchestration
# ─────────────────────────────────────────────────────────────────────────────

def build_outputs_for_variant(
    variant: Variant,
    pbs_len: int = 13,
    max_rtt_len: int = 40,
    require_variant_in_rtt: bool = True,
    codon_table: Optional[dict[str, str]] = None,
) -> list[PegRNAOutput]:
    """Run the full pipeline for one Variant. Returns all (PAM × silent) rows."""
    if codon_table is None:
        codon_table = load_codon_table()

    pams = find_pam_candidates_for_variant(
        variant,
        pbs_len=pbs_len,
        max_rtt_len=max_rtt_len,
        require_variant_in_rtt=require_variant_in_rtt,
    )

    # Compute variant protein change once per variant
    variant_protein = compute_variant_protein_change(variant, codon_table=codon_table)

    rows: list[PegRNAOutput] = []
    for pam in pams:
        # Compute RTT splice context once per PAM (not per silent)
        splice_ctx = _compute_rtt_splice_context(variant, pam)

        # NEW silent_finder API: takes (variant, pam, codon_table)
        silents = find_silent_candidates(
            variant=variant,
            pam=pam,
            codon_table=codon_table,
        )
        for silent in silents:
            rows.append(build_pegrna_output(
                variant, pam, silent,
                variant_protein=variant_protein,
                splice_ctx=splice_ctx,
            ))
    return rows
