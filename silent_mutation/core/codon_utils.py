"""
Codon utilities for silent mutation detection.
 
Pure functions for codon table lookup, translation, synonymy checks,
and frame-aware codon extraction. No dependencies on the io layer —
silent_finder is responsible for unpacking Variant/Transcript objects
and passing primitives in.
 
CDS coordinate convention:
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
 
import csv
 
 
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
 
    with open(csv_path, newline="") as fh:
        reader = csv.DictReader(fh)
        fields = set(reader.fieldnames or [])
        if not {"Codon", "AminoAcid"}.issubset(fields):
            raise ValueError(
                f"Codon table must have 'Codon' and 'AminoAcid' columns. "
                f"Found: {sorted(fields)}"
            )
        table = {
            str(row["Codon"]).upper(): str(row["AminoAcid"]).upper()
            for row in reader
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
