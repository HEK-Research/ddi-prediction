# Results & Notes

Running log of notable findings from the H1/H2/H3 notebooks that aren't obvious from the notebooks alone.

## Data notes

### H1 vs. H2/H3 drug-count discrepancy (2026-09-03)
- H2 and H3 share the exact same sampled pair files (`adverse_atc_sim.parquet` / `adverse_structural_sim.parquet`, etc.): 471,169 adverse pairs / 1,893 unique drugs, 162,788 non-interacting pairs / 1,900 unique drugs.
- H1's pair files (`adverse_biological_overlap_extended.parquet` / `non_interacting_biological_overlap_extended.parquet`) were curated independently by `biological_overlap.ipynb` and differ slightly: 471,286 adverse pairs / 1,894 unique drugs, 163,007 non-interacting pairs / 1,901 unique drugs.
- Root cause: H1 includes 2 drugs (`DB12604`, `DB08884`) absent from H2/H3, while H2/H3 include 1 drug (`DB08840`) absent from H1 -- net +1 drug per population in H1. Not a code bug -- each hypothesis's data-curation notebook filtered/sampled drugs independently (e.g. differing SMILES vs. target-annotation availability upstream).
- Practical impact: when combining embeddings/features across H1/H2/H3 in a downstream model, always filter to the intersection of drug ids embedded by *every* source being combined -- a handful of pairs get dropped due to this mismatch.
