# DDI Prediction Project
Predicting drug-drug interactions using heterogeneous network embeddings and ensemble learning.

### Project Structure
```
ddi-prediction/
|-- notebooks/      # Jupyter notebooks for analysis
|-- notebooks_test/ # In-development, Jupyter notebooks for analysis
|-- src_test/       # In-development, Reusable Python modules 
|-- src/            # Core Reusable Python modules for imports
|-- scripts/        # Standalone scripts
|-- data/           # Data files
    |-- raw/        # Original data (gitignored)
    |-- processed/  # Cleaned data (gitignored)
    |-- README.md   # Data documentation
|-- models/         # Saved models (gitignored)
    |-- README.md   # Model documentation
|-- results/        # Outputs (gitignored)

|-- README.md
|-- rquirements-colab.txt
|-- environment.yml 
|--.gitignore
```

## Setup
### Enrironment Setup
- Conda version: 24.x or higher
- Create environment: "conda env create -f environment.yml"

### Option A: Local Development
```bash
# Clone repository
git clone https://github.com/HEK-Research/ddi-prediction.git
cd ddi-prediction
conda env create -f environment.yml
conda activate ddi-prediction
```

### Option B: Google Colab
1. Open notebook from GitHub in Colab
2. Run setup cell to clone repository
3. Start working! 

## Team Members
- Dr. He (PI)
- Student researchers
    - Ashton Croteau
    - Yesli Linares

## Methods

## License
Research use only