## Examining all drugbank attributes, resulted in the following observations:
1. ### Drug/type: (keep) 15485 small molecule drugs + (discard) 4356biotech drugs.
    

2. ### DDI Lables: 
When drugs are administered concomitantly, will affect its activity or result in adverse effects. 
DDI interactions may be synergistic or antagonistic depending on the physiological effects and mechanism of action of each drug. 
Each <drug> element may have one or more <drug-interaction> elements as children of the <drug-interaction> element. 
        drug-interactions/drug-interaction/drugbank-id 
        drug-interactions/drug-interaction/description

3. ### Only extract these top-level fields: 
    (a) Identity:
        drugbank-id
        name
        cas-number
        unii
        average-mass
        monoisotopic-mass
        groups

    (b) Pharmacological classification 
        classification 
        categories
        atc-codes

    (c) Network edges (further simplify, keep only key identifiers)
        drug-interactions 
        calculated-properties
        experimental-properties
        external-identifiers
        pathways
        reactions
        snp-effects
        snp-adverse-drug-reactions
        targets
        enzymes
        transporters
        carriers

    (d) Pharmacology
        indication
        mechanism-of-action
        pharmacodynamics
        metabolism
        absorption
        half-life
        protein-binding
        route-of-elimination
        volume-of-distribution
        clearance
        toxicity
