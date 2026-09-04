# DDI Feature Completion & Filtering Pipeline

## Overview

This notebook (`weeding_out_unusable_unapproved_drugs.ipynb`) builds the final
drug-drug interaction (DDI) datasets — one **negative** (non-interacting) pair
set and one **adverse-event** (positive interaction) pair set — used for
downstream network-based DDI prediction modeling.

The central design goal:

> Include both **approved** and **unapproved** DrugBank small-molecule drugs
> in the applicability domain, but only keep drugs that have a **complete
> feature set** (SMILES, ATC codes, protein targets). Drugs that still lack
> these core features after extensive multi-source backfilling are dropped.

This produces a middle ground between maximizing drug coverage (applicability
domain) and preserving data quality for model training.

---

## Final Outputs

| File | Description |
|---|---|
| `data/sample/negative_final_df.csv` | Filtered negative DDI pairs, both drugs guaranteed to have complete features |
| `data/sample/adverse_final_df.csv` | Filtered adverse-event DDI pairs, both drugs guaranteed to have complete features |
| `per_drug_df` (in-memory) | Master per-drug feature table used to build both final datasets |

Both final CSVs contain the original pair columns plus per-drug feature
columns suffixed `_1` / `_2` for:

- `smiles`
- `atc_codes_list`
- `target_uniprot_ids`
- `target_fasta_sequences`
- `target_other_ids`
- `target_go_terms`
- `target_pfam_domains`
- `target_pathway_neighbors`

---

## Required vs. Supplemental Features

| Feature | Required for inclusion? | Purpose |
|---|---|---|
| `smiles` | ✅ Required | Chemical structure representation |
| `atc_codes_list` | ✅ Required | Therapeutic classification |
| `target_uniprot_ids` / `target_other_ids` | ✅ Required | Protein target network |
| `target_fasta_sequences` | ⚪ Supplemental | Sequence-based protein features |
| `target_go_terms` | ⚪ Supplemental | Gene Ontology annotation |
| `target_pfam_domains` | ⚪ Supplemental | Protein domain architecture |
| `target_pathway_neighbors` | ⚪ Supplemental | Pathway "guilt-by-association" protein context |

A drug is only discarded at the final step if it's missing SMILES, ATC codes,
**or** targets. The supplemental fields enrich the dataset but don't gate
inclusion.

---

## Pipeline Steps

### 1. Load Raw Pair Data
- Load the raw negative and adverse-event DDI pair CSVs.
- Auto-detect the two DrugBank ID columns in each (`infer_pair_cols`).
- Collect the union of every unique DrugBank ID referenced (`needed_ids`).

### 2. Parse DrugBank XML
- Stream-parse the full DrugBank XML (memory-safe `iterparse`, since the raw
  DB dump is >1 GB and may be zipped).
- For every drug matching a needed ID, extract:
  - SMILES (calculated or experimental properties)
  - ATC codes
  - Target info: UniProt IDs, FASTA sequences, other identifiers, actions
  - Drug groups (approved / unapproved / withdrawn)
- **Filter:** keep only **small-molecule, non-withdrawn** drugs — both
  approved and unapproved are retained at this stage.

### 3. Restrict Pairs to Small-Molecule Drugs
- Build `sm_ids`: the set of small-molecule, non-withdrawn DrugBank IDs.
- Filter both pair datasets so both drugs in each pair belong to `sm_ids`.

### 4. Attach Features & Build Per-Drug Table
- Merge feature columns onto both pair datasets.
- Build `per_drug_df`: one row per unique relevant drug with boolean
  completeness flags — `has_smiles`, `has_atc`, `has_targets`,
  `has_all_features`.
- Report baseline completeness and per-feature gap counts.

### 5. Iterative External Backfilling

A multi-source, multi-pass backfill strategy fills in as many gaps as
possible before any drug is discarded:

1. **PubChem SMILES backfill** — lookup by DrugBank ID xref, fallback to name search.
2. **ChEMBL ATC + target backfill** (parallelized) — resolve ChEMBL molecule ID by
   name, then fetch ATC classifications and mechanism-of-action targets.
3. **Merge round 1** — merge PubChem SMILES + ChEMBL ATC/targets; recompute flags.
4. **UniProt second-pass target search** — for drugs still missing targets,
   search UniProtKB (reviewed only) by drug name or previously scraped gene/target names.
5. **WHO ATC/DDD Index scrape** — best-effort regex-based scrape for drugs still missing ATC codes.
6. **Merge round 2** — merge UniProt targets + WHO ATC codes; recompute flags.
7. **FASTA sequence backfill** (parallelized) — fetch missing sequences for any UniProt accession.
8. **SMPDB/DrugBank pathway-based target inference**
   - Build a shared `PathwayMapper` (loads SMPDB + DrugBank pathway data once).
   - For drugs still missing targets, resolve ChEMBL ID and batch-query
     ChEMBL mechanism/activity data (`min_pchembl=6.0`) for direct targets.
   - Two-hop pathway expansion: target → pathway → co-pathway "neighbor"
     proteins, stored separately as supplemental `target_pathway_neighbors`.
   - Merge newly discovered direct targets into `target_uniprot_ids`; recompute flags.
9. **Batched UniProt sequence/GO/Pfam retrieval** — using `UniprotConverter`,
   batch-fetch sequences, GO terms, and Pfam domains (chunks of 100 IDs) for
   every UniProt accession still missing a sequence, across the whole dataset.
10. **Dataset-wide pathway neighbor expansion** — repeat the two-hop pathway
    expansion across **every** drug's full current target set (not just
    drugs that were missing targets), maximizing pathway-neighbor coverage.

### 6. Final Filtering
- Recompute `has_all_features` (SMILES + ATC + targets required) after all backfilling.
- Discard any drug still missing one of the three core features.
- Filter both pair datasets down to pairs where **both** drugs are in the
  final complete-feature drug set.
- Re-attach the final feature set — including the enrichment columns
  (`target_go_terms`, `target_pfam_domains`, `target_pathway_neighbors`) — to
  both filtered pair datasets via `final_feature_cols`.
- Report approved vs. unapproved drug counts across the complete-feature
  pool, each final dataset, and their union.

### 7. Save Outputs
- Export `negative_final_df` and `adverse_final_df` to CSV in `data/sample/`.

---

## Key Design Decisions

- **Approval status is not a filtering criterion** — only feature completeness is.
  This maximizes the applicability domain while guaranteeing model input quality.
- **Backfilling is exhaustive and layered** — each external source (PubChem,
  ChEMBL, UniProt, WHO ATC, SMPDB pathways) is tried before giving up on a feature,
  minimizing unnecessary drug loss.
- **Pathway neighbors are treated as supplemental, not required** — since they're
  a looser (guilt-by-association) signal rather than direct annotation.
- **All external API calls are parallelized** (`ThreadPoolExecutor`) and batched
  where possible to keep runtime manageable given the large number of unique drugs.

---

## Dependencies

- `pandas`, `numpy`, `requests`
- `xml.etree.ElementTree`, `zipfile`, `pathlib` (standard library)
- Custom module: `src_test/smpdb_protein_pathway.py`
  (`PathwayMapper`, `get_targets_from_chembl_batch`, `UniprotConverter`, `phi_infer_chembl_batch`)

## Required Local Data

- `data/raw/drugbank_full_database_V5.1.14.zip` — full DrugBank XML export
- `data/smpdb_pathways_data_csv/smpdb_proteins.csv.zip` — SMPDB pathway-protein mapping
- Raw negative/adverse pair CSVs (paths configured in Cell 1)