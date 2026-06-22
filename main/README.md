# Silent Mutation Finder — final package

Prime-editing pegRNA tool that finds synonymous (silent) codon substitutions in a
pegRNA's RTT region as MMR-evasion / re-cut-prevention markers. Takes an in-frame
CDS-strand sequence window (WT + Edited), runs PAM finding -> silent finding ->
pegRNA assembly, and optionally ranks pegRNAs by DeepPrime efficiency (Option A:
DeepPrime ranks the therapeutic edit; our silent_finder checks each top pegRNA's RTT).

## Layout
- silent_mutation/core/   types, codon_utils, pam_finder, silent_finder, pegrna_builder
- silent_mutation/io/     sequence_loader (WT/Edited -> Variant), deepprime_runner (Option A),
                          genome_loader (minimal: reverse_complement only)
- silent_mutation/webtool/ server.py (Flask), index.html (UI w/ DeepPrime drill-down)
- data/reference/codon_table.csv

## Run
    cd final
    PYTHONPATH=. python silent_mutation/webtool/server.py    # serves on :8501
DeepPrime ranking requires the DeepPrime-main repo + its deps (tensorflow/torch/ViennaRNA);
without them the standalone path still works (graceful fallback).

## Validated
- EYS c.2528G>A (G843E, 121 candidates) and c.4957dupA (Ser1653 frameshift, 37 candidates),
  extracted from genomic EYS and verified end-to-end.
