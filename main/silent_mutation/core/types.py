"""
Core data structures for silent mutation candidate finding.

Conventions
-----------
- All sequences are stored 5'->3' in CDS strand orientation.
  (Plus-strand genes: same as genome plus strand.
   Minus-strand genes: reverse complement of genome plus strand.)
- All sequence indices into Variant.seq_wt / Variant.seq_ed are 0-based.
- The variant occupies positions [VAR_IDX : VAR_IDX + ref_len] in seq_wt
  and [VAR_IDX : VAR_IDX + alt_len] in seq_ed.
- Window follows SynDesign / DeepPrime convention with FLANK = 60 bp on
  each side of the variant.

Length relationships
--------------------
    len(seq_wt) = FLANK + ref_len + FLANK
    len(seq_ed) = FLANK + alt_len + FLANK

Variant types (by ref_len, alt_len):
    Substitution: ref_len == alt_len  (both 1, 2, or 3)
    Insertion:    ref_len == 0,        alt_len in {1, 2, 3}
    Deletion:     ref_len in {1,2,3},  alt_len == 0
"""

from dataclasses import dataclass, field
from typing import Literal

# SynDesign / DeepPrime standard. Do NOT change without coordinating with
# downstream tools (DeepPrime input format depends on this).

FLANK = 60
VAR_IDX = 60


# ─────────────────────────────────────────────────────────────────────────────
# Variant
# ─────────────────────────────────────────────────────────────────────────────

@dataclass
class Variant:
    """
    A single variant with its surrounding sequence and CDS context.

    Supports substitutions, insertions, and deletions of size 1-3 bp.

    seq_wt is ALWAYS in CDS-strand 5'->3' orientation. For minus-strand genes,
    the io layer takes care of reverse-complementing the genome fetch before
    constructing the Variant. cds_strand is metadata only.
    """

    seq_wt: str
    seq_ed: str
    ref_len: int
    alt_len: int
    cds_strand: Literal["+", "-"]
    cds_frame: int
    cds_start_in_seq: int
    cds_end_in_seq: int
    variant_id: str = ""
    gene_symbol: str = ""
    transcript_id: str = ""
    chrom: str = ""
    genome_pos: int = 0
    ref_allele: str = ""
    alt_allele: str = ""
    codon_lookup: dict = field(default_factory=dict)

    # Splice context — added Step 1 of splice annotation pass
    # near_exon_boundary: True if variant.genome_pos is within ±5 nt of any
    #   CDS exon boundary (potential splice-affecting silent edit risk).
    # nearest_boundary_dist: signed distance (bp) to nearest CDS exon boundary.
    # Both are computed in clinvar_loader when transcript info is available;
    # default to False/0 for variants built without that context.
    near_exon_boundary: bool = False
    nearest_boundary_dist: int = 0

    # ─── derived properties ────────────────────────────────────────────
    @property
    def var_idx(self) -> int:
        return VAR_IDX

    @property
    def variant_type(self) -> Literal["substitution", "insertion", "deletion"]:
        if self.ref_len == 0 and self.alt_len > 0:
            return "insertion"
        if self.ref_len > 0 and self.alt_len == 0:
            return "deletion"
        return "substitution"

    @property
    def is_indel(self) -> bool:
        return self.ref_len != self.alt_len

    @property
    def is_frameshift(self) -> bool:
        """True if the indel changes reading frame (length diff not multiple of 3)."""
        return self.is_indel and ((self.alt_len - self.ref_len) % 3 != 0)

    @property
    def expected_seq_wt_len(self) -> int:
        return FLANK + self.ref_len + FLANK

    @property
    def expected_seq_ed_len(self) -> int:
        return FLANK + self.alt_len + FLANK

    # ─── validation ────────────────────────────────────────────────────
    def __post_init__(self):
        if len(self.seq_wt) != self.expected_seq_wt_len:
            raise ValueError(
                f"seq_wt length {len(self.seq_wt)} != expected "
                f"{self.expected_seq_wt_len} (FLANK={FLANK} + ref_len={self.ref_len} "
                f"+ FLANK={FLANK}). variant_id={self.variant_id!r}"
            )
        if len(self.seq_ed) != self.expected_seq_ed_len:
            raise ValueError(
                f"seq_ed length {len(self.seq_ed)} != expected "
                f"{self.expected_seq_ed_len} (FLANK={FLANK} + alt_len={self.alt_len} "
                f"+ FLANK={FLANK}). variant_id={self.variant_id!r}"
            )

        if self.ref_len not in (0, 1, 2, 3):
            raise ValueError(
                f"ref_len must be 0, 1, 2, or 3 (got {self.ref_len})."
            )
        if self.alt_len not in (0, 1, 2, 3):
            raise ValueError(
                f"alt_len must be 0, 1, 2, or 3 (got {self.alt_len})."
            )
        if self.ref_len == 0 and self.alt_len == 0:
            raise ValueError(
                "At least one of ref_len/alt_len must be > 0."
            )

        if self.cds_strand not in ("+", "-"):
            raise ValueError(f"cds_strand must be '+' or '-' (got {self.cds_strand!r})")
        if self.cds_frame not in (0, 1, 2):
            raise ValueError(f"cds_frame must be 0, 1, or 2 (got {self.cds_frame})")

        if not (0 <= self.cds_start_in_seq <= self.cds_end_in_seq):
            raise ValueError(
                f"Invalid CDS coords: start={self.cds_start_in_seq}, "
                f"end={self.cds_end_in_seq}"
            )
        if self.cds_end_in_seq > self.expected_seq_wt_len:
            raise ValueError(
                f"cds_end_in_seq {self.cds_end_in_seq} exceeds seq_wt length "
                f"{self.expected_seq_wt_len}"
            )

        if self.ref_allele:
            actual_ref = self.seq_wt[VAR_IDX : VAR_IDX + self.ref_len]
            if actual_ref.upper() != self.ref_allele.upper():
                raise ValueError(
                    f"ref_allele mismatch: stored={self.ref_allele!r}, "
                    f"seq_wt slice={actual_ref!r}, variant_id={self.variant_id!r}"
                )
        if self.alt_allele:
            actual_alt = self.seq_ed[VAR_IDX : VAR_IDX + self.alt_len]
            if actual_alt.upper() != self.alt_allele.upper():
                raise ValueError(
                    f"alt_allele mismatch: stored={self.alt_allele!r}, "
                    f"seq_ed slice={actual_alt!r}, variant_id={self.variant_id!r}"
                )

# ─── codon / genome coordinate helpers ─────────────────────────────
    def seq_idx_to_genome_pos(self, i: int) -> int:
        """
        Convert a 0-based seq_wt index to a 1-based genomic position.

        genome_pos is stored as the LEFTMOST plus-strand coordinate anchor
        (VCF POS convention): plus_wt[FLANK] is at genome_pos for every
        variant type (the REF's first base for sub/del; the base immediately
        after the insertion site for insertions).

        The window is FLANK + ref_len + FLANK long. After reverse-complementing
        for a minus-strand gene, plus_wt[FLANK] (at genome_pos) lands at
        seq_wt index FLANK + (ref_len - 1), so the minus-strand mapping needs
        the (ref_len - 1) correction. This term is:
            +1/+2 for 2/3bp variants, 0 for 1bp, and -1 for insertions
                   (ref_len = 0)
        and is exact for ALL of those cases. (An earlier max(ref_len,1)-1
        form was wrong for insertions by one base.)
        """
        if self.cds_strand == "+":
            return self.genome_pos + (i - VAR_IDX)
        return self.genome_pos + (self.ref_len - 1) - (i - VAR_IDX)

    def genome_pos_to_seq_idx(self, gpos: int):
        """
        Inverse of seq_idx_to_genome_pos. Returns the 0-based seq_wt index
        for a given 1-based genomic position, or None if it falls outside
        the seq_wt window.
        """
        if self.cds_strand == "+":
            i = gpos - self.genome_pos + VAR_IDX
        else:
            i = (self.genome_pos - gpos) + VAR_IDX + (self.ref_len - 1)
        if 0 <= i < self.expected_seq_wt_len:
            return i
        return None

    def codon_at_variant(self):
        """
        CodonInfo for the variant's first genomic base, or None if the
        codon_lookup is empty or the variant isn't on a CDS base.
        For indels, this returns the codon at the variant's first ref base
        (or, for insertions where ref_len=0, the codon at the position
        just after the insertion site).
        """
        if not self.codon_lookup:
            return None
        return self.codon_lookup.get(self.genome_pos)

# ─────────────────────────────────────────────────────────────────────────────
# PAM candidate
# ─────────────────────────────────────────────────────────────────────────────

@dataclass
class PamCandidate:
    """
    A single PAM site found near a variant in seq_wt, with derived
    nick / PBS / RTT coordinates.

    All seq_wt-coordinate fields (pam_pos, nick_pos, rtt_start_in_seq,
    rtt_end_in_seq) are 0-based indices into Variant.seq_wt.

    pbs_seq / rtt_seq are stored in PAM-strand 5'->3' orientation,
    matching DeepPrime convention.

      pam_pattern == 'NGG':
        PAM strand == CDS strand
        nick_pos = pam_pos - 3
        pbs/rtt taken directly from seq_wt
      pam_pattern == 'CCN':
        PAM strand != CDS strand
        nick_pos = pam_pos + 6
        pbs/rtt obtained by reverse-complementing the seq_wt slice
    """

    pam_pos: int
    pam_pattern: Literal["NGG", "CCN"]
    nick_pos: int
    pbs_seq: str
    rtt_seq: str
    rtt_start_in_seq: int
    rtt_end_in_seq: int


# ─────────────────────────────────────────────────────────────────────────────
# Silent mutation candidate
# ─────────────────────────────────────────────────────────────────────────────

@dataclass
class SilentCandidate:
    """
    A silent (synonymous) mutation candidate that can be added to a pegRNA
    on top of the intended variant edit.

    Even when the intended variant is an indel, the silent *marker* itself
    is always a substitution within the RTT region.
    """

    pam: PamCandidate
    silent_positions: list[int]
    silent_ref: str
    silent_alt: str
    original_codon: str
    mutated_codon: str
    aa: str
    locale: Literal["PAM", "LHA", "RHA"]
    priority: int


# ─────────────────────────────────────────────────────────────────────────────
# pegRNA output row
# ─────────────────────────────────────────────────────────────────────────────

@dataclass
class PegRNAOutput:
    """
    One pegRNA candidate row, ready to hand off to DeepPrime-coediting.

    Each row = one (variant, PAM, silent_candidate) triple. A single variant
    can produce many rows (multiple PAMs × multiple synonymous alternatives).

    Sequence fields (CDS-strand 5'->3'):
      seq_wt  : the unedited 121bp window (FLANK + ref + FLANK)
      seq_ed  : the variant-edited window (FLANK + alt + FLANK)
      seq_ed_with_silent : seq_ed with the silent substitution applied at
                           silent_positions (still in CDS-strand coords).

    Convention note for downstream:
      - seq_wt and seq_ed match SynDesign's `wtseq` / `edseq` columns and
        feed directly into DeepPrime as Ref_seq / ED_seq.
      - seq_ed_with_silent is the "co-edited" target — what the cell looks
        like when both the intended variant AND the silent marker are
        installed. This is what the pegRNA's RTT actually templates.
    """

    # Identity
    variant_id: str
    gene_symbol: str
    transcript_id: str

    # Variant info
    chrom: str
    genome_pos: int
    ref_allele: str
    alt_allele: str
    cds_strand: Literal["+", "-"]
    variant_type: Literal["substitution", "insertion", "deletion"]
    ref_len: int
    alt_len: int

    # Window sequences (CDS-strand 5'->3')
    seq_wt: str                 # 121bp, variant NOT applied
    seq_ed: str                 # variant applied, silent NOT applied
    seq_ed_with_silent: str     # variant + silent both applied

    # PAM info
    pam_pos: int
    pam_pattern: Literal["NGG", "CCN"]
    nick_pos: int
    pbs_seq: str                # PAM-strand 5'->3'
    rtt_seq_wt: str             # PAM-strand 5'->3', no edits
    rtt_start_in_seq: int
    rtt_end_in_seq: int

    # Silent info
    silent_positions: list[int]   # absolute seq_wt indices, CDS-strand
    silent_ref: str               # ref bases at silent_positions (CDS-strand)
    silent_alt: str               # alt bases at silent_positions (CDS-strand)
    original_codon: str
    mutated_codon: str
    aa: str
    locale: Literal["PAM", "LHA", "RHA"]
    priority: int                 # 1=PAM, 2=LHA, 3=RHA

    # Variant-induced protein change (from codon_lookup at variant.genome_pos)
    # e.g. "D23Y", or "p.D23Y" if you prefer. Empty string for non-CDS or
    # frameshift variants where simple AA substitution doesn't apply.
    variant_protein_change: str = ""
    variant_wt_codon: str = ""
    variant_mut_codon: str = ""
    variant_wt_aa: str = ""
    variant_mut_aa: str = ""
    variant_codon_index: int = 0

    # Variant-level splice context (copied from Variant for per-row analysis)
    variant_near_exon_boundary: bool = False
    variant_nearest_boundary_dist: int = 0

    # RTT-level splice context
    # cds_bases_in_rtt: how many of the RTT's bases are CDS (in codon_lookup).
    #   Lower means the RTT spills into intron / UTR → fewer silent options.
    # rtt_covers_intron: True if RTT contains any non-CDS base (intron/UTR).
    # splice_junction_codons_in_rtt: count of codons touched by RTT whose
    #   three bases are not contiguous in the genome (splice-junction codons).
    cds_bases_in_rtt: int = 0
    rtt_covers_intron: bool = False
    splice_junction_codons_in_rtt: int = 0
