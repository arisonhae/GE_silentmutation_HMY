"""
silent_mutation.io.genome_loader  (minimal shim for the webtool deliverable)

The full io/genome_loader.py pulls in pyfaidx + pyranges for genome-FASTA work,
which the sequence-input demo never touches and which are NOT installed in the
lab 'dprime' env. The webtool compute path (pam_finder, deepprime_runner) only
needs reverse_complement, so this lightweight version keeps `final/` importable
with zero heavy genomics dependencies.

Behaviour is identical to the full module's reverse_complement.
"""

_COMPLEMENT = str.maketrans("ACGTNacgtn", "TGCANtgcan")


def reverse_complement(seq: str) -> str:
    return seq.translate(_COMPLEMENT)[::-1]
