## Drug Representation via ATC Code Path encoding
Official website: https://atcddd.fhi.no/atc_ddd_index/

### Background
Each drug must be represented as a fixed-length numerical vector suitable for machine learning. 
This document desdribes how ATC codes are used to construct a multi-hot path encoding for each drug. 

For example, 
| Metformin (A10BA02) |
|---|---|---|
| Level 1 | one letter | A |
| Level 2 | two digits | A10 |
| Level 3 | one letter | A10B |
| Level 4 | one letter | A10BA |

The codes do not carry standalone meaning and are not shared across each upper levels. 

### ATC-based Hierarchical path similarity
#### 1. Rationale
A one-hot encoding at a single ATC level loses the hierarchical relationship between drugs. 
Two drugs in the same chemical subclass (Level 4) are more pharmacoclogically similar than two drugs that merely share the same anatomical group (Level 1). 

#### 2. Levels Used
Levels 2, 3, and 4 are used for encoding. 
Level 1 is excluded because it provides only 14 coarse categories with limited discriminative power. 
Level 5 (individual drug substance) is excluded because it would encode drug identify rather than pharmacological class membership, which is the signal of interest for DDI prediction. 

Level 1 Anatomical main group   14
Level 2 Therapeutic subgroup    ~94
Level 3 Pharmacological subgroup    ~271
Level 4 Chemical Subgroup   ~939

**Example**
Given ATC code, each drug maps to a set of ATC level-lists:
|drug name | ATC code | ATC Level list|
|"metformin"| (A10BA02)|  ["A10", "A10B", "A10BA"]|
|"glipizide"| (A10BB07) | ["A10", "A10B", "A10BB"]|

#### 3. Encoding Procedure
For each drug, the encoding proceeds as follows: 


#### 5. Handling Multiple ATC Codes


#### 6. Similarity Metrics for ATC-based Hierarchical Path Encoding
Given two ATC codes with level representations, depth weight is defined. 
Match indicator: a path match at depth d only counts if all shallower levels also matched 


**Assigning Wegihts to Bit Positions**
For the weighted Tanimoto to work correctly, each position in the drug encoding vector must know its weight. 
This is determined by the ATC code stored at that index position: the code length maps directly to ATC levek, which in tern maps to weight: 

ATC_CODE_LENGTH ATC_Level   Example Weight
3               Level 2     A10     1
4               Level 3     A10B    2
5               Level 4     A10BA   3

Since the index is an ordered list of ATC codes (index_codes), the weight vector is built once from that list and reused for every similarity computation. Each position i in the weight vector corresponds to the same position i in the drug encoding vector. 

#### 7. Limitations
* The encoding captures pharmacological class membership but not drug-specific properties such as molecular structure, physicochemical parameters, or binding affinity. 
* Drugs without ATC code cannot be encoded with this method and require an alternative representation. 
* The WHO ATC index is updated annually, the current encoding is based on the 2026 edition. Future updates may add new codes, requiring re-indexing.
* The encoding is sparse and high-dimensional relative to the information content. 

#### 8. Related Files
File                                            Description
src_test/scrape_atc.py                          Scrapes ATC classes (Level 1-4) from atcddd.fhi.no and saves to CSV.
data/atc_code_output/WHO_ATC_codes_<date>.csv   Complete list of ATC codes with level labels
data/atc_code_output/rds_cache/                 Per-root pickle cache files. Not tracked in Git (.gitignore)



