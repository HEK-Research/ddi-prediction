# LLM Classification Runs - System Architecture & Batch Processing

## Overview

This directory implements a production-grade batch processing system for classifying Drug-Drug Interactions (DDIs) using Google's Gemini LLM. The system processes large datasets efficiently through concurrent batch operations, stores each run independently with comprehensive outputs, and maintains strict traceability through timestamped run folders.

## File System Architecture

```
lmm_classifcation_runs/
├── README.md                              # User guide for using the classification system
├── SYSTEM_ARCHITECTURE.md                 # This file - technical architecture
├── run_20260709_202349/                   # Example run folder (timestamp-based naming)
│   ├── 01_prompt_documentation.txt        # System prompt & configuration used
│   ├── 02_classification_results.csv      # Main results: all 1500 DDI classifications
│   ├── 03_statistics_summary.json         # Structured stats for analysis
│   ├── 04_category_distribution.json      # Detailed breakdown by category
│   ├── 05_boundary_validation_analysis.txt# Quality assurance notes
│   ├── 06_final_run_summary.txt           # Human-readable summary
│   ├── 07_execution_metadata.json         # Machine-readable metadata
│   ├── 08_distribution_bar_chart.png      # Visualization: bar chart
│   └── 09_distribution_pie_chart.png      # Visualization: pie chart
└── run_YYYYMMDD_HHMMSS/                   # Future runs follow same structure
    ├── [9 standardized output files]
    └── [all self-contained within run folder]
```

## Batch Processing Strategy

### Sample Size: 1500 DDIs per Run

Each classification run processes **1500 randomly sampled Drug-Drug Interactions** from the full DrugBank dataset. This size balances:
- **Computational efficiency**: Reasonable execution time with API costs
- **Statistical significance**: Large enough for reliable category distribution analysis
- **Manageability**: Fits comfortably in memory and API quotas

### Batch Size: 50 Concurrent Requests

The 1500 samples are processed in **30 batches of 50 concurrent requests each**:

```
Total samples: 1500 DDIs
Batches: 30 (1500 ÷ 50 = 30)
Concurrency per batch: 50 simultaneous API calls
Processing flow:

Batch 1: Items 0-49     [50 concurrent requests] → Wait 0.5-1.0s delay
Batch 2: Items 50-99    [50 concurrent requests] → Wait 0.5-1.0s delay
Batch 3: Items 100-149  [50 concurrent requests] → Wait 0.5-1.0s delay
...
Batch 30: Items 1450-1499 [50 concurrent requests] → Done
```

### Rate Limiting & Retry Logic

The system implements intelligent rate limiting:

1. **Inter-batch delays**: 0.5-1.0 second delay between batches prevents quota exhaustion
2. **Exponential backoff retries**: 
   - Rate limit errors (429): 5s, 10s, 20s, 40s, 80s waits
   - Other errors: 1s, 2s, 4s, 8s, 16s waits
3. **Error handling**: Failed classifications recorded with error details in CSV

### Async/Concurrent Architecture

```python
# Pseudocode flow:
for each batch of 50:
    tasks = [classify_interaction_async(desc) for desc in batch_descriptions]
    results = await asyncio.gather(*tasks)  # All 50 run concurrently
    save_batch_results()
    await sleep(inter_batch_delay)  # Rate limit protection
```

**Benefits**:
- 50 DDIs classified simultaneously instead of sequentially
- Reduces total execution time from ~25 minutes (serial) to ~3-5 minutes (parallel)
- API quota efficiently utilized without overwhelming the service

## Run Organization & Naming

### Timestamped Run Folders

Each execution creates a **unique timestamped folder** using the format:
```
run_YYYYMMDD_HHMMSS
```

**Example**: `run_20260709_202349` = July 9, 2026 at 20:23:49 (8:23:49 PM)

**Benefits**:
- **Immutability**: Each run is self-contained and never overwritten
- **Traceability**: Exact execution time embedded in folder name
- **Sortability**: Lexicographic ordering = chronological ordering
- **Uniqueness**: Collision-proof with second-level precision

### Run Independence

Each run folder is **completely independent** and contains:
- Complete input configuration (system prompt, model settings)
- Full results (all 1500 classifications)
- Comprehensive statistics and metadata
- Visualization outputs

This enables:
- Running multiple classification jobs without conflicts
- Comparing results across different LLM configurations
- Archiving historical runs for reproducibility
- Batch processing different sample sets simultaneously

## 9 Standardized Output Files

Every run generates exactly 9 output files, numbered `01-09`:

| # | File | Type | Purpose |
|---|------|------|---------|
| 1 | `01_prompt_documentation.txt` | Text | System prompt & configuration used |
| 2 | `02_classification_results.csv` | CSV | All 1500 DDI classifications with reasoning |
| 3 | `03_statistics_summary.json` | JSON | Structured statistics (machine-readable) |
| 4 | `04_category_distribution.json` | JSON | Detailed category breakdown with samples |
| 5 | `05_boundary_validation_analysis.txt` | Text | Quality assurance & edge case analysis |
| 6 | `06_final_run_summary.txt` | Text | Human-readable summary report |
| 7 | `07_execution_metadata.json` | JSON | Machine-readable run metadata |
| 8 | `08_distribution_bar_chart.png` | PNG | Bar chart visualization (1000×600px) |
| 9 | `09_distribution_pie_chart.png` | PNG | Pie chart visualization (900×700px) |

**Standardization Benefits**:
- Consistent file naming across all runs
- Predictable file locations for automation
- Easy comparison between runs
- Clear organization for both humans and scripts

## Classification Configuration

### LLM Model
- **Model**: `gemini-2.5-flash`
- **Temperature**: 0.0 (deterministic for reproducibility)
- **Response format**: Structured JSON with Pydantic validation

### 15-Category System

**13 Adverse Effect Categories**:
1. Bleeding/Hemorrhage Risk
2. Dosage Adjustment Required
3. Organ Toxicity
4. CNS/Neurological Effects
5. Metabolic/Absorption Interference
6. Cardiovascular Effects
7. Hematological Effects
8. Fluid & Electrolyte Imbalance
9. Clinical Metabolic Dysregulation
10. Immunological/Hypersensitivity
11. Therapeutic Efficacy Alteration
12. Gastrointestinal
13. General/Unspecified Adverse Effects

**2 Beneficial/Neutral Categories**:
14. Positive/Synergistic
15. No Significant Interaction

### System Prompt

Each run documents the exact system prompt used (see `01_prompt_documentation.txt`), which includes:
- Category definitions with clinical explanations
- Critical boundary rules (e.g., "Bleeding ≠ Hematological")
- Priority hierarchy for multi-effect interactions
- Example-driven classification guidance

## Data Flow Diagram

```
Raw DrugBank Data
      ↓
   [Sample 1500 random DDIs]
      ↓
   [Split into 30 batches of 50]
      ↓
   ┌─────────────────────────────────┐
   │  Batch 1: Async Process (50)    │
   │  - 50 concurrent API calls      │
   │  - Each: "Classify this DDI"    │
   │  - LLM returns category + reason │
   │  - Pydantic validation          │
   └─────────────────────────────────┘
   ↓ [Wait 0.5-1.0s]
   ┌─────────────────────────────────┐
   │  Batch 2: Async Process (50)    │
   │  ...                             │
   └─────────────────────────────────┘
   ↓ [Repeat 30 times]
      ↓
   [Collect all 1500 results]
      ↓
   ┌─────────────────────────────────┐
   │  Generate 9 Output Files:       │
   │  - CSVs, JSONs, TXTs, PNGs      │
   │  - Save to run_YYYYMMDD_HHMMSS/ │
   └─────────────────────────────────┘
```

## Usage Examples

### Load Results from a Run
```python
import pandas as pd
import json

run_id = "run_20260709_202349"

# Load classifications
results = pd.read_csv(f"{run_id}/02_classification_results.csv")
print(results.head())

# Load statistics
with open(f"{run_id}/03_statistics_summary.json") as f:
    stats = json.load(f)
    print(f"Total samples: {stats['total_samples']}")
    print(f"Distribution: {stats['category_distribution']}")
```

### Compare Multiple Runs
```bash
# List all runs chronologically
ls -d run_* | sort

# Compare statistics across runs
for run in run_*; do
    echo "=== $run ==="
    cat $run/06_final_run_summary.txt
done
```

## API Cost & Performance Metrics

### Execution Time
- **Total runtime**: ~3-5 minutes per 1500-sample run
- **Per-DDI classification**: ~150-200ms average (including retries)
- **Batch overhead**: Minimal with concurrent processing

### API Costs (Gemini 2.5 Flash)
- **Input tokens**: ~50,000 per run (30 batches × system prompt re-use)
- **Output tokens**: ~7,500 per run (1500 classifications × 5 tokens)
- **Total cost**: Minimal on paid tier

### Resource Usage
- **Memory**: ~50-100 MB per batch (1500 samples in memory)
- **Network**: ~2-5 MB (mostly JSON structures)
- **Disk**: ~10-20 MB per run folder

## Quality Assurance

Each run includes quality checks:

1. **Pydantic validation**: Every classification validated against schema
2. **Boundary verification**: Critical categories audited (see `05_boundary_validation_analysis.txt`)
3. **Error tracking**: Failed items recorded with error details
4. **Statistics verification**: Distribution summaries computed and cross-checked
5. **Metadata logging**: All run parameters documented for reproducibility

## Future Runs & Scaling

The architecture supports:
- **Multiple concurrent runs**: Each uses separate folder, no conflicts
- **Different sample sizes**: Easily configurable (1500 is default)
- **Custom batching**: Adjust batch size for different API quotas
- **Historical analysis**: All past runs preserved with metadata
- **Comparative studies**: Run different LLM models/prompts, compare results

## File Navigation Quick Reference

```
To analyze results:
  → Read 02_classification_results.csv (primary data)
  → View 06_final_run_summary.txt (executive summary)
  → Check 03_statistics_summary.json (metrics)
  → View PNG charts for presentations

To understand the run:
  → Read 01_prompt_documentation.txt (what LLM was told)
  → Check 07_execution_metadata.json (run details)

To debug or verify:
  → Review 05_boundary_validation_analysis.txt (QA notes)
  → Check 04_category_distribution.json (detailed breakdown)
```

---

**Last Updated**: 2026-07-09  
**System Version**: 15-Category Classification v1.0  
**Architecture**: Async batch processing with rate limiting and error recovery
