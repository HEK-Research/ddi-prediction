## Drug Representation via ATC Code Path encoding
Official website: https://atcddd.fhi.no/atc_ddd_index/

### Background
Each drug must be represented as a fixed-length numerical vector suitable for machine learning. 
This document desdribes how ATC codes are used to construct a multi-hot path encoding for each drug. 

For example, Metformin (A10BA02) 
Level 1-one letter: A
Level 2-two digits: A10
Level 3-one letter: A10B
Level 4-one letter: A10BA

The codes do not carry standalone meaning and are not shared across each upper levels. 

### Multi-Hot Path Encoding
#### 1. Rationale
A one-hot encoding at a single ATC level loses the hierarchical relationship between drugs. 
Two drugs in the same chemical subclass (Level 4) are more pharmacoclogically similar than two drugs that merely share the same anatomical group (Level 1). 
Multi-hot path encoding addresses this by encoding the full taxonomic path of a drug, capturing similarity at every level simultaneously. 

#### 2. Levels Used
Levels 2, 3, and 4 are used for encoding. 
Level 1is excluded because it provides only 14 coarse categories with limited discriminative power. 
Level 5 (individual drug substance) is excluded because it would encode drug identify rather than pharmacological class membership, which is the signal of interest for DDI prediction. 

The fixed index in constructed from the complete WHO ATC reference (scraped from atcddd.fhi.no), not from the drug dataset itself, to ensure consistent vector dimensions across any dataset: 

Level 1 Anatomical main group   14
Level 2 Therapeutic subgroup    ~94
Level 3 Pharmacological subgroup    ~271
Level 4 Chemical Subgroup   ~939
Total vector dimension = 94 + 271 + 939 = 1304

#### 3. Encoding Procedure
For each drug, the encoding proceeds as follows: 
* Initialize a zero vector of length 1,304.
* Retrieve all ATC codes assigned to the drug (a drug may have multiple codes for different indications or routes of administration). 
* For each ATC code, extract its level 2, 3, and 4 ancestors. 
* Set the corresponding positions in the vector to 1 (logical OR across all codes).

The result is a binary vector in which a 1 at position i means the drug belongs to ATC category i via at least one of its assgiend codes. 

A consistent encoding across all datasets requires a fixed mapping from ATC code to vector index, built once from the complete WHO ATC reference. The index is constructed as follows: 
* Scrape all Level 1-4 ATC codes from atcddd.fhi.no using scrape_atc.py
* Filter to Levels 2, 3, and 4 only
* Sort codes alphabetically within each level (L2 first, then L3, then L4)
* Assign a sequential integer index starting from 0. 

This index is saved as a CSV file and loaded at encoding time. 
Any drug from any dataset - including drugs not present in the training set - can be encoded into this same fixed-dimensional space. 

**Note**: If a drug's ATC code contains an ancestor not present in the WHO 2026 index, that ancestor is silently ignored. This is rare and does not affect the majority of DrugBank entires. 

Property |              Description
Fixed dimensionality    All drugs map to a vector of length 1,304 regardless of how many ATC codes they carry
Binary                  Each element is 0 or 1, no magnitude information is encoded
Sparse                  Most drugs belong to a small fraction of the 1,304 categories
Hierarchy-aware         Drugs sharing a Level 4 ancestor will share at least 3 bits

#### 4. Example
Two example drugs illustrates the procedure. For clarity, only a subset of relevant index positions is shown:

ATC_CODE|    Level|   Drug1(A10BA02)|  Drug2(A10BB07)|  Drug3(C10AA05)|  Drug4(C10AB02)
A10         2       1               1               0               0
C10         2       0               0               1               1
A10B        3       1               1               0               0
C10A        3       0               0               1               1
A10BA       4       1               0               0               0
A10BB       4       0               1               0               0
C10AA       4       0               0               1               0
C10AB       4       0               0               0               1

#### 5. Handling Multiple ATC Codes
**The OR Operation**
OR is a logical operation that returns 1 if at least one input is 1, and 0 only if both inputs are 0:

Using addition instead of OR would assign a value of 2 to any bit shared by two coes, breaking the binary nature of the vector and artificially over-weighting shared ancestors in any similarity computation. 

* Case 1: Non-overlapping ATC codes (different Level 1 branches)
When a drug has codes from completely differnt Level 1 branches, there is no shared ancestry - the merged vector is simply the union of all bits, with no overlap to resolve. This accurately reflects that the drugs has dual therapeutic roles spanning two completely different pharmacological systems. 

* Case 2: Partially overlapping ATC codes (shared ancestors)
When two codes share ancestors at higher levels but diverge at a lower level, OR correctly duplicates the shared bits. 

**Implementation**
    vec = np.zeros(1304, dtype=int)
    for atc_code in drug_atc_codes:
        for ancestor in get_ancestors(atc_code): # return L2, L3, L4
            vec[atc_to_index[ancestor]] = 1 
            

#### 6. Similarity Metrics for ATC Path Encoding
The Tanimoto coefficient is identical to the Jaccard Similarity for binary vectors, and is the standard metric in cheminformatics for comparing binary fingerprints. It is the recommended metric for ATC path encodings. 

Its key advantages is that it ignores shared zeros, whcih is essential for sparse binary data. 

The standard Tanimoto treats all bits equally. However, a shared Level 4 category is pharmacologically more informative than a shared level 2 category. A weighted variant assigns greater weight to deeper levels: 

Level 2 weight = 1, Level 3 weight = 2, Level 4 weight = 3

The weighted variant penalizes shallow similarity more heavily and rewards deep Level 4 overlap more strongly, better reflecting true pharmacological relatedness. 

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



