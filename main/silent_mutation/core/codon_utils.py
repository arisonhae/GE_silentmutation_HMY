"""
Codon utilities for silent mutation detection.
 
Pure functions for codon table lookup, translation, synonymy checks,
and frame-aware codon extraction. No dependencies on the io layer —
silent_finder is responsible for unpacking Variant/Transcript objects
and passing primitives in.
 
CDS coordinate convention (established in io layer v4):
    - seq is stored 5'→3' in CDS strand orientation
      (minus-strand transcripts are already RC-normalized by genome_loader)
    - cds_start, cds_end are 0-based half-open [start, end) within seq
    - cds_frame ∈ {0, 1, 2} = number of nucleotides preceding the first
      complete codon in the CDS region. The first complete codon begins at
      seq[cds_start + cds_frame].
"""
 
from __future__ import annotations
 
from functools import lru_cache
from pathlib import Path
from typing import Optional
 
import pandas as pd
 
 
# Default location: <project_root>/data/reference/codon_table.csv
_DEFAULT_CODON_TABLE_PATH = (
    Path(__file__).resolve().parents[2] / "data" / "reference" / "codon_table.csv"
)
 
_VALID_BASES = frozenset("ACGT")
 
 
# ---------------------------------------------------------------------------
# Codon table loading
# ---------------------------------------------------------------------------
 
@lru_cache(maxsize=4)
def load_codon_table(path: Optional[str] = None) -> dict[str, str]:
    """
    Load the codon → amino acid mapping from CSV.
 
    Expected CSV format (header required):
        Codon,AminoAcid
        ATA,I
        ...
        TAA,*
 
    Stop codons are represented as '*'. Result is cached per path.
 
    Args:
        path: Optional path to codon table CSV. Defaults to
            data/reference/codon_table.csv relative to project root.
 
    Returns:
        Dict mapping uppercase 3-letter codon strings to single-letter AA codes.
 
    Raises:
        FileNotFoundError: If the CSV file does not exist.
        ValueError: If the table is malformed or incomplete (≠64 codons).
    """
    csv_path = Path(path) if path is not None else _DEFAULT_CODON_TABLE_PATH
 
    if not csv_path.exists():
        raise FileNotFoundError(f"Codon table not found: {csv_path}")
 
    df = pd.read_csv(csv_path)
 
    if not {"Codon", "AminoAcid"}.issubset(df.columns):
        raise ValueError(
            f"Codon table must have 'Codon' and 'AminoAcid' columns. "
            f"Found: {list(df.columns)}"
        )
 
    table = {
        str(row["Codon"]).upper(): str(row["AminoAcid"]).upper()
        for _, row in df.iterrows()
    }
 
    if len(table) != 64:
        raise ValueError(
            f"Codon table must contain exactly 64 codons, got {len(table)}"
        )
 
    # Sanity check: every key must be a valid 3-letter ACGT codon
    for codon in table:
        if len(codon) != 3 or not set(codon).issubset(_VALID_BASES):
            raise ValueError(f"Invalid codon in table: {codon!r}")
 
    return table
 
 
# ---------------------------------------------------------------------------
# Translation
# ---------------------------------------------------------------------------
 
def translate_codon(codon: str, table: Optional[dict[str, str]] = None) -> str:
    """
    Translate a 3-letter codon to a single-letter amino acid code.
 
    Args:
        codon: 3-letter DNA codon (case-insensitive). Must be ACGT only.
        table: Optional pre-loaded codon table. If None, loads default.
 
    Returns:
        Single-letter amino acid code. Stop codons return '*'.
 
    Raises:
        ValueError: If codon length ≠ 3 or contains non-ACGT characters.
    """
    if len(codon) != 3:
        raise ValueError(f"Codon must be 3 nucleotides, got {len(codon)}: {codon!r}")
 
    codon_upper = codon.upper()
 
    if not set(codon_upper).issubset(_VALID_BASES):
        raise ValueError(
            f"Codon contains non-ACGT characters: {codon!r}"
        )
 
    if table is None:
        table = load_codon_table()
 
    return table[codon_upper]
 
 
def is_synonymous(
    codon1: str,
    codon2: str,
    table: Optional[dict[str, str]] = None,
) -> bool:
    """
    Check whether two codons encode the same amino acid.
 
    Note: This is a "dumb" comparison — stop↔stop returns True. The caller
    (e.g. silent_finder) is responsible for any policy around excluding
    stop-codon-involved synonymy from silent marker candidates.
 
    Args:
        codon1: First 3-letter DNA codon.
        codon2: Second 3-letter DNA codon.
        table: Optional pre-loaded codon table. If None, loads default.
 
    Returns:
        True if both codons translate to the same amino acid (or both stop).
 
    Raises:
        ValueError: If either codon is malformed.
    """
    if table is None:
        table = load_codon_table()
 
    return translate_codon(codon1, table) == translate_codon(codon2, table)
 
 
# ---------------------------------------------------------------------------
# Frame-aware codon extraction
# ---------------------------------------------------------------------------
 
def get_codon_at_position(
    seq: str,
    position: int,
    cds_start: int,
    cds_frame: int,
) -> Optional[tuple[str, int, int]]:
    """
    Extract the codon containing a given nucleotide position.
 
    Frame logic:
        - cds_frame is the offset of the first complete codon within the CDS.
        - The first complete codon begins at seq[cds_start + cds_frame].
        - For a position p inside the CDS, its offset from the first codon
          start is: rel = p - (cds_start + cds_frame)
          The position-in-codon is rel % 3, and the codon starts at p - (rel % 3).
 
    Positions before cds_start + cds_frame (i.e. in the partial leading
    codon, if any) are treated as out-of-frame and return None — silent
    mutation analysis only operates on complete in-frame codons.
 
    Positions outside [cds_start, cds_end) are also out of scope and return
    None. Note: this function does not know cds_end; the caller must
    upstream-filter, OR rely on the codon falling off the end of seq
    (which raises a length check below).
 
    Args:
        seq: Sequence in CDS strand orientation (uppercase recommended).
        position: 0-based absolute position within seq.
        cds_start: 0-based start of CDS within seq (inclusive).
        cds_frame: Frame offset {0, 1, 2} of first complete codon in CDS.
 
    Returns:
        (codon_str, codon_start_in_seq, position_in_codon) or None if the
        requested position falls outside an extractable in-frame codon.
        position_in_codon ∈ {0, 1, 2} indicates which slot of the codon
        the requested position occupies.
 
    Raises:
        ValueError: If cds_frame is not in {0, 1, 2}.
    """
    if cds_frame not in (0, 1, 2):
        raise ValueError(f"cds_frame must be 0, 1, or 2, got {cds_frame}")
 
    # Position must be at or after the first complete codon start
    first_codon_start = cds_start + cds_frame
    if position < first_codon_start:
        return None
 
    # Position must be inside seq
    if position >= len(seq):
        return None
 
    rel = position - first_codon_start
    pos_in_codon = rel % 3
    codon_start = position - pos_in_codon
 
    # Bounds check: codon must fit within seq
    if codon_start + 3 > len(seq):
        return None
 
    codon = seq[codon_start : codon_start + 3].upper()
 
    # Defensive: if the codon contains non-ACGT (e.g. N from masked region),
    # treat as un-extractable rather than raising — biologically expected at
    # assembly gaps, and silent_finder should skip these candidates.
    if not set(codon).issubset(_VALID_BASES):
        return None
 
    return codon, codon_start, pos_in_codon