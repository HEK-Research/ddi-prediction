# LLM Classification Runs

This directory stores organized outputs from LLM-based DDI classification runs using the 15-category system.

## Directory Structure

Each classification run creates a timestamped subdirectory: `run_YYYYMMDD_HHMMSS/`

### Run Output Files

Each run folder contains the following standardized outputs:

#### 1. `01_prompt_documentation.txt`
- **Purpose**: Complete documentation of the system prompt used
- **Contents**:
  - System configuration (model, temperature, parameters)
  - Full 15-category system prompt with all definitions
  - Category structure explanation
  - Key improvements over previous systems
  - Output artifacts description

#### 2. `02_classification_results.csv`
- **Purpose**: Primary results file with all classifications
- **Contents**: CSV with columns:
  - `drug1_name`: First drug in interaction
  - `drug2_name`: Second drug in interaction
  - `description`: DDI description from DrugBank
  - `classification`: Assigned category (one of 15)
  - `Reason`: LLM's reasoning for classification

#### 3. `03_statistics_summary.json`
- **Purpose**: Structured statistics for programmatic use
- **Contents**:
  - Run ID and timestamp
  - Model and temperature settings
  - Total samples
  - Number of categories (15)
  - Category distribution counts
  - Grouping statistics (adverse vs. beneficial)

#### 4. `04_category_distribution.json`
- **Purpose**: Detailed breakdown by category with sample examples
- **Contents**:
  - For each of 15 categories:
    - Count and percentage
    - 5 representative sample interactions
    - Drug pairs with descriptions

#### 5. `05_boundary_validation_analysis.txt`
- **Purpose**: LLM verification of classification accuracy and boundary cases
- **Contents**:
  - Sample interactions from key categories
  - LLM analysis of boundary accuracy
  - Potential misclassification identification
  - Validation against category definitions

#### 6. `06_final_run_summary.txt`
- **Purpose**: Human-readable summary of the complete run
- **Contents**:
  - Run ID, timestamp, and model info
  - Final statistics with full category breakdown
  - New categories vs. original categories comparison
  - File index and usage notes

#### 7. `07_execution_metadata.json`
- **Purpose**: Machine-readable metadata for run tracking
- **Contents**:
  - Run ID and timestamps
  - Model configuration
  - Sample counts and category breakdowns
  - Output directory path
  - List of all generated files

#### 8. `08_distribution_bar_chart.png`
- **Purpose**: Bar chart visualization
- **Format**: PNG image (1000x600px)
- **Contents**:
  - Count of DDIs in each of the 15 categories
  - High-quality image suitable for presentations and reports

#### 9. `09_distribution_pie_chart.png`
- **Purpose**: Pie chart visualization
- **Format**: PNG image (900x700px)
- **Contents**:
  - Percentage distribution across all 15 categories
  - High-quality image suitable for presentations and reports

## 15-Category Classification System

### Adverse Effects (13 categories)
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

### Beneficial/Neutral (2 categories)
14. Positive/Synergistic
15. No Significant Interaction

## Key Features

- **Deterministic**: Temperature set to 0.0 for consistent classifications
- **Structured Output**: Validated via Pydantic schema
- **Well-Documented**: Complete prompt and reasoning captured
- **Rate-Limited**: Implements exponential backoff for API stability
- **Organized**: Timestamp-based folders for easy tracking
- **Reproducible**: All configuration and prompts documented
- **Visualized**: PNG plots for presentations and reports
- **Comprehensive**: 9 output artifacts covering all aspects of the run

## Usage Examples

### Load Classification Results
```python
import pandas as pd
df = pd.read_csv('run_20260709_143022/02_classification_results.csv')
```

### Access Run Metadata
```python
import json
with open('run_20260709_143022/07_execution_metadata.json') as f:
    metadata = json.load(f)
    print(f"Total samples: {metadata['sample_size']}")
    print(f"Categories: {metadata['num_categories']}")
```

### View Statistics
```python
import json
with open('run_20260709_143022/03_statistics_summary.json') as f:
    stats = json.load(f)
    print(json.dumps(stats, indent=2))
```

### View Plots
The PNG image files can be viewed directly:
- `run_20260709_143022/08_distribution_bar_chart.png` - Bar chart
- `run_20260709_143022/09_distribution_pie_chart.png` - Pie chart

Open in any image viewer or embed in documents. High-quality PNG format suitable for presentations and publications.

## Analysis Workflow

1. **Run the notebook** (cells execute in order)
2. **Outputs are generated** automatically with timestamps
3. **Browse the run folder** to review results
4. **Compare across runs** using metadata files
5. **Integrate downstream** using CSV results or JSON metadata

## Notes

- Each run is independent and timestamped
- Results are final (temperature=0.0 means no randomness)
- All artifacts needed for reproducibility are included
- CSV files can be directly imported to tools like Excel, pandas, or databases
- JSON files maintain structured format for programmatic access
- PNG images are high-quality and suitable for presentations and publications
