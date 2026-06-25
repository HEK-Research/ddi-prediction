# DDI Prediction Project
Predicting drug-drug interactions using heterogeneous network embeddings and ensemble learning.

## Development-Branch (Stage I, Chemical Similarities)
Main task: 
1. Gather all drugbank compounds that met the filtering standard 
2. Using random sampling, create 10 (n=100) subsets for method testing
3. Develop chem_similarities.py, that process each drug and calculate pair-wise similarity scores
4. Develop chem_network.py, that reads the similarity matrix and construct network
5. Develop chem_embedding.py, that applies node2vec algorithm to extract drug's embedding representation
6. Develop chem_rf_model.py, that use network embedding as input features to train a binary predictive model 



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