"""
smpdb_protein_pathway.py
================================================================================

WHAT THIS SCRIPT DOES
---------------------
This script builds "pathway-level protein sets" for drugs, which are used to
detect drug-drug interactions (DDIs) that operate through shared biology rather
than shared direct targets.

The core insight: two drugs can interact not because they bind the same protein,
but because they both perturb the same biological pathway — even through entirely
different proteins. This script provides the tools to detect that kind of
indirect overlap.


KEY CONCEPTS
------------
Pathway
    A named group of proteins that cooperate in a biological process,
    e.g. "Blood Coagulation" or "Vitamin K Metabolism". A protein can
    belong to many pathways at once.

Prot(s)  [get_proteins_for_pathway]
    The atomic lookup: given a pathway ID, return every protein in it.
    All other operations are built on top of this.

Π(u)  [get_pathways_for_protein]
    The inverse lookup: given a protein, return every pathway it belongs to.
    This is the first "hop" in the two-hop expansion.

Φ_native(d)  [phi_native]
    The native pathway profile of a drug.
    Starts from the drug's curated SMPDB pathway annotations (from DrugBank)
    and expands each pathway to its full protein membership. The result is
    the complete protein landscape of the biological space the drug was
    *annotated* to act in. Think of it as: "populate every room the drug
    is known to occupy with all the proteins in those rooms".

Φ_infer(d; l)  [phi_infer]
    The inferred pathway profile of a drug.
    Starts from the drug's known protein interactions (targets, enzymes,
    transporters — the L_l set) and does a two-hop expansion:
        Hop 1: for each interacting protein, find all its pathways
        Hop 2: for each of those pathways, collect all member proteins
    The result is every protein biologically reachable from the drug's
    targets through shared pathway membership. Think of it as: "starting
    from who the drug touches, find everyone else in the same rooms".

DDI Signal
    Two drugs are at risk of interacting when their expanded protein sets
    overlap. Overlap in phi_native means they act in the same biological
    context. Overlap in phi_infer means their targets share pathway
    neighborhoods — a more mechanistic signal.


DATA SOURCES
------------
All four are queried automatically and their results are merged per ID:

  DrugBank XML   Local file (~800 MB). Best source for SMPDB pathway data.
                 Scanned once per call using streaming parse to stay memory-
                 efficient. Pass xml_path= to any function to enable it.

  SMPDB REST     smpdb.ca/api/v1 — currently returns HTTP 403 for automated
                 requests; DrugBank is the reliable fallback for SMP* IDs.

  Reactome REST  reactome.org/ContentService — no credentials needed.
                 The most reliable live source; covers R-HSA-* pathway IDs.

  ChEMBL REST    ebi.ac.uk/chembl — no credentials needed.
                 Provides protein→pathway cross-references and drug→target
                 bioactivity data (MoA + experimental assay records).


QUICK START
-----------
# 1. Get all proteins in a pathway
from smpdb_protein_pathway import get_proteins_for_pathway
proteins = get_proteins_for_pathway("R-HSA-159740")
# → ['P00734', 'P00740', 'P00742', ...]

# 2. Get all pathways a protein belongs to
from smpdb_protein_pathway import get_pathways_for_protein
pathways = get_pathways_for_protein("P00734")
# → ['R-HSA-140837', 'R-HSA-159740', ...]

# 3. Φ_native — expand a drug's curated SMPDB pathway list to proteins
from smpdb_protein_pathway import phi_native
native_proteins = phi_native(["SMP0000278", "SMP0000765"])
# → {'P00451', 'P00734', ...}

# 4. Φ_infer — two-hop expansion from a drug's target proteins
from smpdb_protein_pathway import phi_infer
result = phi_infer(["P00734", "P00742"])
result["pathways"]    # all pathways those targets collectively sit in
result["neighbors"]   # all proteins that share any of those pathways

# 5. ChEMBL pipeline — drug ID directly to inferred pathway neighborhood
from smpdb_protein_pathway import phi_infer_chembl
result = phi_infer_chembl("CHEMBL1536", min_pchembl=6)  # warfarin
result["proteins"]    # UniProt targets found in ChEMBL
result["pathways"]    # pathways those targets sit in
result["neighbors"]   # full expanded protein neighborhood

# 6. Annotate proteins with sequence, GO terms, Pfam domains
from smpdb_protein_pathway import get_sequence, get_go_terms, get_pfams
get_sequence("P00734")    # → 'MAHVRGLQLPG...'
get_go_terms("P00734")    # → ['GO:0004252', 'GO:0005509', ...]
get_pfams("P00734")       # → ['PF00089', 'PF09396']

# 7. Batch annotation for many proteins at once
from smpdb_protein_pathway import UniprotConverter
conv = UniprotConverter(["P00734", "P00742", "P00533"])
conv.uniprot_to_sequence()    # → {'P00734': 'MAHVR...', ...}
conv.uniprot_to_GO_terms()    # → {'P00734': ['GO:...'], ...}
conv.uniprot_to_pfams()       # → {'P00734': ['PF...'], ...}

# 8. Efficient batch lookups for many IDs at once
from smpdb_protein_pathway import PathwayMapper
mapper = PathwayMapper(xml_path="drugbank_full_database.xml")  # optional
mapper.get_proteins_by_pathway(["SMP0000278", "R-HSA-159740"])
mapper.get_pathways_by_protein(["P00734", "P00533"])


NOTES
-----
- All UniProt accessions follow the standard format: P00734, Q9Y6K9, etc.
- Pathway ID prefixes: SMP* = SMPDB, R-HSA-* = Reactome, hsa* = KEGG,
  WP* = WikiPathways
- min_pchembl controls ChEMBL bioactivity stringency:
    None  → all records with any measured potency
    5     → IC50/Ki ≤ 10 µM  (broad)
    6     → IC50/Ki ≤  1 µM  (drug-like, recommended)
    7     → IC50/Ki ≤ 100 nM (potent/selective)
"""

import gzip
import xml.etree.ElementTree as ET
import zipfile
from concurrent.futures import ThreadPoolExecutor
import pandas as pd
import requests


# Every tag in the DrugBank XML carries this namespace prefix, e.g. <ns:drug>, <ns:name>.
# We must include it when searching, otherwise ET.find/findall return nothing.
NS = "{http://www.drugbank.ca}"

# Base URLs so they're easy to swap out
_SMPDB_BASE    = "https://www.smpdb.ca/api/v1"
_REACTOME_BASE = "https://reactome.org/ContentService"
_CHEMBL_BASE   = "https://www.ebi.ac.uk/chembl/api/data"
_UNIPROT_BASE  = "https://rest.uniprot.org"


# ── Prot(s), Π(u), Φ_native, Φ_infer — the four core framework operations ──
# These are the simple entry points. Use them when you just need a quick answer
# without spinning up the full PathwayMapper.


def get_proteins_for_pathway(pathway_id: str, xml_path: str | None = None) -> list[str]:
    """Prot(s) — the atomic lookup: given a pathway, return all proteins in it.

    Conceptually, a pathway is a named biological process (e.g. "Blood
    Coagulation") with a defined set of participating proteins. This function
    answers: "Who is in this room?"  It queries all available sources and
    merges the results, so you get the most complete answer regardless of
    which database originally defined the pathway.

    This is the primitive operation — phi_native and phi_infer both call it
    internally as part of their expansions.

    Args:
        pathway_id: Any supported pathway identifier.
                    SMP* → SMPDB,  R-HSA-* → Reactome,  hsa* → KEGG,
                    WP* → WikiPathways.
        xml_path:   Path to the DrugBank full-database XML. Leave as None to skip.

    Returns:
        Sorted list of UniProt accessions. Empty list if nothing is found.

    Example:
        proteins = get_proteins_for_pathway("SMP0000278")
        # ['P00451', 'P00488', 'P00734', ...]
    """
    return PathwayMapper(xml_path).get_proteins_by_pathway(pathway_id)[pathway_id]


def get_pathways_for_protein(protein_id: str, xml_path: str | None = None) -> list[str]:
    """Π(u) — the inverse lookup: given a protein, return all pathways it belongs to.

    Conceptually, a single protein can participate in many pathways — thrombin
    (P00734) appears in blood coagulation, fibrinolysis, platelet activation, and
    more. This function answers: "Which rooms does this protein appear in?"

    This is the first hop used by phi_infer: you start from a drug's targets,
    call this on each one to find all their pathways, then expand those pathways
    outward to discover the full biological neighborhood.

    Note: this function only does the first hop. To get the full two-hop expansion
    (proteins → pathways → neighbor proteins), use phi_infer instead.

    Args:
        protein_id: A UniProt accession, e.g. "P00734".
        xml_path:   Path to the DrugBank full-database XML. Leave as None to skip.

    Returns:
        Sorted list of pathway IDs spanning all sources (SMP*, R-HSA-*, hsa*, WP*).

    Example:
        pathways = get_pathways_for_protein("P00734")
        # ['R-HSA-140837', 'R-HSA-159740', 'SMP0000278', ...]
    """
    return PathwayMapper(xml_path).get_pathways_by_protein(protein_id)[protein_id]


def phi_native(pathway_ids: str | list[str], xml_path: str | None = None) -> set[str]:
    """Φ_native(d) — the native pathway profile of a drug.

    Conceptual idea:
        DrugBank curates a list of SMPDB pathways for each drug — the biological
        processes the drug is known to participate in. This function takes that
        list and expands it: for each pathway, collect every protein that
        participates in it, then union everything together.

        The result is the complete protein landscape of the biological space
        the drug was *annotated* to act in. Think of it as:
            "Given the rooms this drug is known to occupy,
             who else is in those rooms?"

    Formal definition:
        Φ_native(d) = ⋃_{s ∈ P_native(d)} Prot(s)

        where P_native(d) is the drug's curated SMPDB pathway list from DrugBank.

    DDI use:
        If Drug A and Drug B both expand to overlapping phi_native sets, they
        share the same biological context — a candidate DDI signal via pathway
        co-membership rather than direct target overlap.

    Args:
        pathway_ids: The drug's curated SMPDB pathway IDs from DrugBank, e.g.
                     ["SMP0000278", "SMP0000765"]. You get these by parsing
                     the <pathways> block of the drug's DrugBank record.
        xml_path:    Path to the DrugBank full-database XML. Strongly recommended
                     for SMP* IDs since SMPDB's live API currently blocks access.

    Returns:
        Set of UniProt accessions — every protein found across all the input
        pathways, merged into one flat set.

    Example:
        profile = phi_native(["SMP0000278", "SMP0000765"])
        # {'P00451', 'P00488', 'P00734', 'P00740', ...}
    """
    mapper  = PathwayMapper(xml_path)
    ids     = mapper._normalize(pathway_ids)
    results = mapper.get_proteins_by_pathway(ids)
    return {p for proteins in results.values() for p in proteins}


def phi_infer(
    protein_ids: str | list[str],
    xml_path: str | None = None,
) -> dict:
    """Φ_infer(d; l) — the inferred pathway profile of a drug's protein interaction set.

    Conceptual idea:
        Rather than starting from annotated pathways (like phi_native does),
        this function starts from the proteins the drug *actually touches* —
        its targets, enzymes, transporters — and asks: given these proteins,
        which other proteins end up biologically reachable through shared
        pathway membership?

        Think of it as:
            "Starting from who the drug touches, find everyone else in
             the same biological rooms."

        This is a two-hop graph traversal:
            Hop 1 (proteins → pathways):  for each input protein, find every
                pathway it participates in. This is Π(u).
            Hop 2 (pathways → proteins):  for each of those pathways, collect
                all their member proteins. This is Prot(s).
            The union of all collected proteins is the Φ_infer neighborhood.

    Formal definition:
        Φ_infer(d; l) = ⋃_{u ∈ L_l(d)}  ⋃_{s ∈ Π(u)}  Prot(s)

        where L_l(d) is the drug's protein interaction set at layer l:
            L0 = direct pharmacological targets only
            L1 = targets + metabolic enzymes + transporters + carriers
            L_ChEMBL = targets resolved from ChEMBL MoA + bioactivity records
        This function is layer-agnostic — pass whichever set is appropriate.

    DDI use:
        If Drug A's phi_infer neighborhood overlaps with Drug B's direct targets
        (or vice versa), then Drug A's biological reach extends into Drug B's
        mechanism — a strong mechanistic DDI signal.

    Args:
        protein_ids: A UniProt accession or list of accessions — the drug's
                     L_l protein interaction set, e.g. ["P00734", "P08183"].
        xml_path:    Path to the DrugBank full-database XML. Leave as None to skip.

    Returns:
        A dict with two keys:
          "pathways"  — sorted list of all pathway IDs collected in Hop 1
          "neighbors" — sorted list of UniProt accessions collected in Hop 2
                        (the full inferred neighborhood; includes input proteins)

    Example:
        result = phi_infer(["P00734", "P08183"])
        print(result["pathways"][:3])
        # ['R-HSA-140837', 'R-HSA-159740', 'SMP0000278']
        print(len(result["neighbors"]))
        # 214
    """
    mapper   = PathwayMapper(xml_path)
    ids      = mapper._normalize(protein_ids)

    # Step 1: for each protein u, call Π(u) to find every pathway it belongs to.
    pathway_map = mapper.get_pathways_by_protein(ids)
    all_pathways: set[str] = set()
    for pw_list in pathway_map.values():
        all_pathways.update(pw_list)

    if not all_pathways:
        return {"pathways": [], "neighbors": []}

    # Step 2: for each pathway s found above, call Prot(s) to get its members,
    # then take the union across all pathways — that’s the Φ_infer neighborhood.
    protein_map   = mapper.get_proteins_by_pathway(sorted(all_pathways))
    all_neighbors: set[str] = set()
    for members in protein_map.values():
        all_neighbors.update(members)

    return {
        "pathways":  sorted(all_pathways),
        "neighbors": sorted(all_neighbors),
    }


# ── ChEMBL target expansion ───────────────────────────────────────────────────
# These functions let you go beyond DrugBank by pulling every target ChEMBL
# has curated for a drug, then feeding that expanded protein set into phi_infer.


def _chembl_targets_to_uniprot(target_ids: list[str]) -> list[str]:
    """Resolve a list of ChEMBL target IDs to UniProt accessions.

    Uses the ChEMBL target __in batch filter, chunked to stay under URL limits.
    Only single-protein targets contribute an accession; complex/organism-level
    targets are skipped automatically because their components have no accession.

    Args:
        target_ids: ChEMBL target IDs, e.g. ["CHEMBL3243", "CHEMBL4523"].

    Returns:
        Sorted, deduplicated list of UniProt accessions.
    """
    uniprot: set[str] = set()
    chunk_size = 50
    for i in range(0, len(target_ids), chunk_size):
        chunk   = target_ids[i : i + chunk_size]
        targets = PathwayMapper._chembl_paginate(
            f"{_CHEMBL_BASE}/target"
            f"?target_chembl_id__in={','.join(chunk)}&format=json&limit=100"
        )
        for target in targets:
            for component in target.get("target_components", []):
                acc = component.get("accession", "").strip()
                if acc:
                    uniprot.add(acc)
    return sorted(uniprot)


def get_targets_from_chembl(
    molecule_chembl_id: str,
    min_pchembl: float | None = None,
) -> list[str]:
    """Return UniProt accessions for all curated targets of a ChEMBL molecule.

    Pulls from two sources and merges them:

    1. **Mechanism of action targets** — ChEMBL's curated MoA table, high-confidence,
       always included regardless of min_pchembl.
    2. **Bioactivity targets** — all assay records where a potency value was measured
       (pChEMBL value present). Optionally filtered by a minimum pChEMBL threshold.

    pChEMBL is -log10(activity in molar), so:
        pChEMBL ≥ 5  →  IC50/Ki ≤ 10 µM   (broad sweep)
        pChEMBL ≥ 6  →  IC50/Ki ≤  1 µM   (typical drug-like)
        pChEMBL ≥ 7  →  IC50/Ki ≤ 100 nM  (potent)

    Args:
        molecule_chembl_id: ChEMBL molecule ID, e.g. "CHEMBL1536" (warfarin).
        min_pchembl:        If set, only include bioactivity records where
                            pChEMBL >= this value. None = include all records
                            that have any measured pChEMBL value.

    Returns:
        Sorted list of UniProt accessions. Empty list if none are found.

    Example:
        proteins = get_targets_from_chembl("CHEMBL1536", min_pchembl=6)
        # ['P00734', 'P01325', ...]
    """
    target_ids: set[str] = set()

    # 1 — mechanism of action targets: ChEMBL's manually curated primary targets.
    #     These are always included regardless of any potency threshold.
    mechs = PathwayMapper._chembl_paginate(
        f"{_CHEMBL_BASE}/mechanism"
        f"?molecule_chembl_id={molecule_chembl_id}&format=json&limit=200"
    )
    for m in mechs:
        tid = m.get("target_chembl_id", "").strip()
        if tid:
            target_ids.add(tid)

    # 2 — bioactivity targets: every single-protein assay record that has a
    #     measured pChEMBL value (pChEMBL = -log10(activity in molar)).
    #     Using pchembl_value__isnull=false avoids records with no numeric readout.
    activities = PathwayMapper._chembl_paginate(
        f"{_CHEMBL_BASE}/activity"
        f"?molecule_chembl_id={molecule_chembl_id}"
        f"&target_type=SINGLE+PROTEIN"
        f"&pchembl_value__isnull=false"
        f"&format=json&limit=1000"
    )
    for a in activities:
        if min_pchembl is not None:
            pval = a.get("pchembl_value")
            try:
                # Skip records that fall below the requested potency cutoff.
                if pval is None or float(pval) < min_pchembl:
                    continue
            except (TypeError, ValueError):
                continue  # malformed pchembl_value field — skip safely
        tid = a.get("target_chembl_id", "").strip()
        if tid:
            target_ids.add(tid)

    if not target_ids:
        return []

    return _chembl_targets_to_uniprot(sorted(target_ids))


def phi_infer_chembl(
    molecule_chembl_id: str,
    min_pchembl: float | None = None,
    xml_path: str | None = None,
) -> dict:
    """Full pipeline: ChEMBL drug ID → protein set → inferred pathway neighborhood.

    Conceptual idea:
        DrugBank's target lists are carefully curated but conservative — they
        only include well-established, textbook-level interactions. ChEMBL
        captures the broader experimental record: thousands of published assays
        testing the drug against many proteins, including off-targets and secondary
        pharmacology that never made it into DrugBank.

        This function uses ChEMBL as the protein source instead of DrugBank, then
        feeds those proteins straight into phi_infer to get the full two-hop
        pathway expansion. It answers:
            "Starting from every protein this drug has been tested against
             (above a potency threshold), what is its full pathway neighborhood?"

        It assembles the protein set from two ChEMBL layers:
            1. MoA targets — ChEMBL's manually curated mechanism-of-action table.
               High confidence. Always included regardless of min_pchembl.
            2. Bioactivity targets — every single-protein assay record that has
               a measured pChEMBL value. Filtered by min_pchembl if set.

    Formal definition:
        L_ChEMBL(d) = MoA_targets(d) ∪ { u : pChEMBL(d,u) ≥ threshold }
        Φ_infer(d; ChEMBL) = phi_infer( L_ChEMBL(d) )

    pChEMBL guide (= -log10 of the IC50/Ki in molar):
        None → all records with any numeric potency measurement
        5    → ≤ 10 µM   broad sweep, catches weak binders
        6    → ≤  1 µM   drug-like potency  ← recommended starting point
        7    → ≤ 100 nM  potent and selective interactions only

    Args:
        molecule_chembl_id: ChEMBL molecule ID, e.g. "CHEMBL1536" (warfarin).
        min_pchembl:        Potency cutoff applied to bioactivity records.
                            MoA targets are always included regardless of this value.
        xml_path:           Path to DrugBank XML for SMPDB pathway coverage.
                            Leave as None to use Reactome + ChEMBL sources only.

    Returns:
        Dict with three keys:
          "proteins"  — UniProt accessions resolved from ChEMBL (the L_ChEMBL set)
          "pathways"  — pathway IDs those proteins collectively sit in (Hop 1)
          "neighbors" — UniProt accessions that share any of those pathways (Hop 2)

    Example:
        result = phi_infer_chembl("CHEMBL1536", min_pchembl=6)
        print(result["proteins"])
        # ['P11473', ...]  — proteins ChEMBL found for warfarin
        print(result["pathways"][:3])
        # ['R-HSA-196791', 'R-HSA-383280', ...]
        print(len(result["neighbors"]))
        # 70
    """
    proteins = get_targets_from_chembl(molecule_chembl_id, min_pchembl=min_pchembl)
    if not proteins:
        print(f"No ChEMBL targets found for {molecule_chembl_id}")
        return {"proteins": [], "pathways": [], "neighbors": []}

    inferred = phi_infer(proteins, xml_path=xml_path)
    return {
        "proteins":  proteins,
        "pathways":  inferred["pathways"],
        "neighbors": inferred["neighbors"],
    }


# ── UniProt annotation helpers ────────────────────────────────────────────────
# One accession in, one result out. For batching many accessions at once,
# use UniprotConverter below.


def get_sequence(uniprot_id: str) -> str | None:
    """Fetch the canonical amino-acid sequence for a UniProt accession.

    Uses the FASTA endpoint — no ID-mapping job needed since we already have
    a canonical accession.

    Args:
        uniprot_id: A UniProt accession, e.g. "P00734".

    Returns:
        The sequence as a plain string (no header or line breaks), or None if
        the accession couldn't be fetched.

    Example:
        seq = get_sequence("P00734")
        print(seq[:20])
        # 'MAHVRGLQLPGCLALAALCP'
    """
    resp = requests.get(f"{_UNIPROT_BASE}/uniprotkb/{uniprot_id}.fasta")
    if resp.status_code != 200:
        print(f"Couldn't fetch sequence for {uniprot_id} (HTTP {resp.status_code})")
        return None

    # FASTA format: first line is the header (>sp|...), everything else is sequence
    lines = resp.text.splitlines()
    return "".join(lines[1:])


def _fetch_uniprot_xrefs(uniprot_id: str, fields: str) -> list[dict]:
    """Fetch cross-reference entries for a single accession from the UniProt entry API.

    Uses GET /uniprotkb/{acc}?fields=...&format=json — no async job needed since
    we already have canonical UniProt accessions.

    Args:
        uniprot_id: A UniProt accession, e.g. "P00734".
        fields:     Comma-separated UniProt fields, e.g. "accession,go,xref_pfam".

    Returns:
        List of cross-reference dicts from uniProtKBCrossReferences, or [] on error.
    """
    resp = requests.get(
        f"{_UNIPROT_BASE}/uniprotkb/{uniprot_id}",
        params={"fields": fields, "format": "json"},
    )
    if resp.status_code != 200:
        print(f"UniProt entry lookup failed for {uniprot_id} (HTTP {resp.status_code})")
        return []
    return resp.json().get("uniProtKBCrossReferences", [])


def get_go_terms(uniprot_id: str) -> list[str]:
    """Return GO term IDs annotated on a UniProt accession.

    Biological process, molecular function, and cellular component terms are all
    included — filter by GO prefix yourself if you only want one category.

    Args:
        uniprot_id: A UniProt accession, e.g. "P00734".

    Returns:
        List of GO term IDs, e.g. ["GO:0004252", "GO:0005509", ...].
        Empty list if the accession isn't found or has no GO annotations.

    Example:
        terms = get_go_terms("P00734")
        print(terms[:3])
        # ['GO:0001525', 'GO:0002576', 'GO:0004252']
    """
    xrefs = _fetch_uniprot_xrefs(uniprot_id, fields="accession,go")
    return [x["id"] for x in xrefs if x.get("database") == "GO"]


def get_pfams(uniprot_id: str) -> list[str]:
    """Return Pfam domain IDs annotated on a UniProt accession.

    Args:
        uniprot_id: A UniProt accession, e.g. "P00734".

    Returns:
        List of Pfam domain IDs, e.g. ["PF00089", ...].
        Empty list if no Pfam entries are annotated.

    Example:
        domains = get_pfams("P00734")
        print(domains)
        # ['PF00089', 'PF09396']
    """
    xrefs = _fetch_uniprot_xrefs(uniprot_id, fields="accession,xref_pfam")
    return [x["id"] for x in xrefs if x.get("database") == "Pfam"]


# ── UniprotConverter — batch UniProt annotation fetching ─────────────────────

class UniprotConverter:
    """Fetch sequence, GO terms, and Pfam domains for a batch of UniProt accessions.

    When you have a large set of proteins — for example, the neighbors returned
    by phi_infer — you often want to characterise them all at once. This class
    fetches all three annotation types in parallel using a thread pool, so
    querying 100 proteins takes roughly the same wall-clock time as querying 5.

    The three annotation types:
        Sequence   — the canonical amino acid sequence (from UniProt FASTA).
                     Used for structural analysis or sequence similarity.
        GO terms   — Gene Ontology annotations describing what the protein does:
                     molecular function, biological process, cellular component.
                     Useful for functional enrichment analysis on the neighbor set.
        Pfam IDs   — Protein family domain annotations. Useful for identifying
                     which protein families appear in a pathway neighborhood and
                     whether two drugs affect similar domain architectures.

    For a single protein, use the standalone helpers instead:
        get_sequence(), get_go_terms(), get_pfams()

    Args:
        uniprot_id: A single accession or list of accessions, e.g. "P00734"
                    or ["P00734", "P00533", "P11712"].

    Example:
        conv = UniprotConverter(["P00734", "P00533"])
        seqs = conv.uniprot_to_sequence()   # {'P00734': 'MAHVR...', ...}
        go   = conv.uniprot_to_GO_terms()   # {'P00734': ['GO:0004252', ...], ...}
        pfam = conv.uniprot_to_pfams()      # {'P00734': ['PF00089', ...], ...}
    """

    def __init__(self, uniprot_id: str | list[str]):
        self.uniprot_ids: list[str] = (
            [uniprot_id] if isinstance(uniprot_id, str) else list(uniprot_id)
        )

    def _batch_xrefs(self, field: str, db_name: str) -> dict[str, list[str]]:
        """Fetch cross-references of one database type for all accessions in parallel."""
        def fetch(acc: str) -> tuple[str, list[str]]:
            xrefs = _fetch_uniprot_xrefs(acc, fields=f"accession,{field}")
            return acc, [x["id"] for x in xrefs if x.get("database") == db_name]

        workers = max(1, min(len(self.uniprot_ids), 20))
        with ThreadPoolExecutor(max_workers=workers) as pool:
            return dict(pool.map(fetch, self.uniprot_ids))

    def uniprot_to_sequence(self) -> dict[str, str]:
        """Return the canonical sequence for each accession, fetched in parallel.

        Returns:
            Dict of {accession: sequence_string}. Accessions that fail are omitted.
        """
        def fetch(acc: str) -> tuple[str, str | None]:
            return acc, get_sequence(acc)

        workers = max(1, min(len(self.uniprot_ids), 20))
        with ThreadPoolExecutor(max_workers=workers) as pool:
            return {acc: seq for acc, seq in pool.map(fetch, self.uniprot_ids) if seq is not None}

    def uniprot_to_GO_terms(self) -> dict[str, list[str]]:
        """Return GO term IDs for each accession in the batch.

        Returns:
            Dict of {accession: [go_term_ids]}.
        """
        return self._batch_xrefs("go", "GO")

    def uniprot_to_pfams(self) -> dict[str, list[str]]:
        """Return Pfam domain IDs for each accession in the batch.

        Returns:
            Dict of {accession: [pfam_ids]}.
        """
        return self._batch_xrefs("xref_pfam", "Pfam")


# ── PathwayMapper — efficient batch pathway lookups across four databases ─────

class PathwayMapper:
    """Efficient batch mapper: proteins ↔ pathways, across all four databases.

    The standalone functions (get_proteins_for_pathway, get_pathways_for_protein,
    phi_native, phi_infer) are convenient for one-off queries. When you need to
    process many IDs at once — e.g. mapping hundreds of proteins from a full drug
    target list — this class is significantly more efficient because it:

        - Scans the DrugBank XML *once* for the entire batch (not once per ID)
        - Fires all REST requests to SMPDB, Reactome, and ChEMBL in parallel
        - Merges results from all four sources into a single answer per ID

    The two public methods are mirrors of each other:
        get_proteins_by_pathway(ids)  — given pathway IDs, return their proteins
        get_pathways_by_protein(ids)  — given protein accessions, return their pathways

    Both accept a single ID or a list, and both return a dict keyed by the input IDs.
    Mixed pathway namespaces (SMP*, R-HSA-*, hsa*, WP*) are handled in one call.

    Args:
        xml_path: Path to the DrugBank full-database XML (.xml, .gz, or .zip).
                  Strongly recommended for SMP* pathway coverage since SMPDB's
                  live API currently blocks automated access. Pass None to skip.

    Example:
        mapper = PathwayMapper("drugbank_full_database.xml")

        # pathway → proteins
        result = mapper.get_proteins_by_pathway(["SMP0000278", "R-HSA-159740"])
        # {"SMP0000278": ["P00734", ...], "R-HSA-159740": ["P00734", ...]}

        # protein → pathways
        result = mapper.get_pathways_by_protein(["P00734", "P00533"])
        # {"P00734": ["R-HSA-140837", "SMP0000278", ...], "P00533": [...]}
    """

    _CHEMBL_PATHWAY_DBS = {"Reactome", "KEGG", "WikiPathways", "PathWhiz", "PharmGKB_Pathway"}
    _CHEMBL_CHUNK = 50   # max IDs per ChEMBL __in filter to stay under URL length limits
    _MAX_WORKERS  = 20   # thread pool size for parallelised REST calls

    def __init__(self, xml_path: str | None = None):
        self.xml_path = xml_path

    # ── Utilities ─────────────────────────────────────────────────────────────

    @staticmethod
    def _normalize(ids: str | list[str]) -> list[str]:
        """Coerce a single string or a list to a deduplicated list (order preserved)."""
        return [ids] if isinstance(ids, str) else list(dict.fromkeys(ids))

    @staticmethod
    def _open_xml(xml_path: str):
        """Return a readable file handle for a DrugBank XML file.

        Handles the three formats DrugBank distributes: plain .xml, .gz, and .zip.
        """
        if xml_path.endswith(".zip"):
            zf        = zipfile.ZipFile(xml_path, "r")
            xml_files = [f for f in zf.namelist() if f.endswith(".xml")]
            if not xml_files:
                raise FileNotFoundError("No XML file found inside the zip archive.")
            # Open the first XML member directly from the zip without extracting it.
            return zf.open(xml_files[0])
        if xml_path.endswith(".gz"):
            return gzip.open(xml_path, "rb")
        # Plain uncompressed XML.
        return open(xml_path, "rb")

    @staticmethod
    def _is_smpdb_id(pid: str) -> bool:
        return pid.upper().startswith("SMP")

    @staticmethod
    def _is_reactome_id(pid: str) -> bool:
        return pid.upper().startswith("R-")

    @staticmethod
    def _parallel_get(
        ids: list[str],
        url_fn,
        parse_fn,
        n_workers: int,
        label: str = "",
    ) -> dict[str, set[str]]:
        """Fire GET requests for all IDs in parallel and collect results.

        Args:
            ids:       IDs to query — each becomes one request.
            url_fn:    Called with an ID, returns the URL string.
            parse_fn:  Called with the response JSON, returns a set of strings.
            n_workers: Thread pool size.
            label:     Source name used in error messages (e.g. "SMPDB").

        Returns:
            {id: set_of_results}. 404s produce an empty set silently.
        """
        result: dict[str, set[str]] = {pid: set() for pid in ids}

        def fetch(pid: str) -> tuple[str, set[str]]:
            resp = requests.get(url_fn(pid), headers={"Accept": "application/json"})
            # 400 and 404 both mean "not found for this ID format" — return empty silently.
            # ChEMBL in particular returns 400 (not 404) for unsupported ID formats.
            if resp.status_code in (400, 404):
                return pid, set()
            if resp.status_code != 200:
                if label:
                    print(f"{label} lookup failed for {pid} (HTTP {resp.status_code})")
                return pid, set()
            return pid, parse_fn(resp.json())

        # pool.map preserves order and returns (pid, result) tuples;
        # dict.update accepts an iterable of (key, value) pairs directly.
        with ThreadPoolExecutor(max_workers=n_workers) as pool:
            result.update(pool.map(fetch, ids))
        return result

    @staticmethod
    def _chembl_paginate(url: str) -> list[dict]:
        """Walk through all pages of a ChEMBL paginated response and return the full list."""
        results  = []
        next_url: str | None = url
        while next_url:
            resp = requests.get(next_url)
            if resp.status_code != 200:
                break
            data = resp.json()
            # ChEMBL responses always have exactly two top-level keys:
            # "page_meta" (pagination info) and one payload key (e.g. "activities",
            # "mechanisms", "targets"). We grab whichever key isn’t page_meta.
            payload_key = next(k for k in data if k != "page_meta")
            results.extend(data[payload_key])
            # page_meta.next is a ready-to-use URL for the next page, or null when done.
            next_url = (data.get("page_meta") or {}).get("next")
        return results

    # ── DrugBank (single XML scan for the whole batch) ────────────────────────

    def _drugbank_pathways_by_proteins(self, protein_ids: set[str]) -> dict[str, set[str]]:
        """Scan the DrugBank XML once and collect SMPDB pathway IDs for every requested protein."""
        if not self.xml_path:
            return {pid: set() for pid in protein_ids}

        result: dict[str, set[str]] = {pid: set() for pid in protein_ids}
        with self._open_xml(self.xml_path) as fh:
            for _event, elem in ET.iterparse(fh, events=["end"]):
                if elem.tag != f"{NS}drug":
                    continue
                for pathway in elem.findall(f"{NS}pathways/{NS}pathway"):
                    smpdb_id = (pathway.findtext(f"{NS}smpdb-id") or "").strip()
                    if not smpdb_id:
                        continue
                    for uid in pathway.findall(f"{NS}enzymes/{NS}uniprot-id"):
                        acc = (uid.text or "").strip()
                        if acc in protein_ids:
                            result[acc].add(smpdb_id)
                # Release the parsed element immediately to keep memory flat.
                # Without this the entire ~800 MB XML would accumulate in RAM.
                elem.clear()
        return result

    def _drugbank_proteins_by_pathways(self, pathway_ids: set[str]) -> dict[str, set[str]]:
        """Scan the DrugBank XML once and collect UniProt IDs for every requested pathway."""
        if not self.xml_path:
            return {pid: set() for pid in pathway_ids}

        result: dict[str, set[str]] = {pid: set() for pid in pathway_ids}
        with self._open_xml(self.xml_path) as fh:
            for _event, elem in ET.iterparse(fh, events=["end"]):
                if elem.tag != f"{NS}drug":
                    continue
                for pathway in elem.findall(f"{NS}pathways/{NS}pathway"):
                    smpdb_id = (pathway.findtext(f"{NS}smpdb-id") or "").strip()
                    if smpdb_id not in pathway_ids:
                        continue
                    for uid in pathway.findall(f"{NS}enzymes/{NS}uniprot-id"):
                        acc = (uid.text or "").strip()
                        if acc:
                            result[smpdb_id].add(acc)
                # Same streaming memory trick — clear after each <drug> block.
                elem.clear()
        return result

    # ── SMPDB (parallelised — no batch endpoint) ──────────────────────────────

    def _smpdb_pathways_by_proteins(self, protein_ids: list[str]) -> dict[str, set[str]]:
        """Look up SMPDB pathways for all proteins in parallel."""
        w = max(1, min(len(protein_ids), self._MAX_WORKERS))
        return self._parallel_get(
            protein_ids,
            url_fn   = lambda pid: f"{_SMPDB_BASE}/pathways/uniprot/{pid}",
            parse_fn = lambda j: {e.get("smpdb_id", "").strip() for e in j if e.get("smpdb_id", "").strip()},
            n_workers = w,
            label    = "SMPDB pathway",
        )

    def _smpdb_proteins_by_pathways(self, pathway_ids: list[str]) -> dict[str, set[str]]:
        """Look up SMPDB proteins for all SMP* pathway IDs in parallel."""
        smpdb_ids = [pid for pid in pathway_ids if self._is_smpdb_id(pid)]
        result: dict[str, set[str]] = {pid: set() for pid in pathway_ids}
        if smpdb_ids:
            w = max(1, min(len(smpdb_ids), self._MAX_WORKERS))
            result.update(self._parallel_get(
                smpdb_ids,
                url_fn   = lambda pid: f"{_SMPDB_BASE}/pathways/{pid}/proteins",
                parse_fn = lambda j: {e.get("uniprot_id", "").strip() for e in j if e.get("uniprot_id", "").strip()},
                n_workers = w,
                label    = "SMPDB protein",
            ))
        return result

    # ── Reactome (parallel for both directions) ───────────────────────────────

    def _reactome_pathways_by_proteins(self, protein_ids: list[str]) -> dict[str, set[str]]:
        """Look up Reactome pathways for all proteins in parallel via the entity endpoint."""
        w = max(1, min(len(protein_ids), self._MAX_WORKERS))
        return self._parallel_get(
            protein_ids,
            url_fn   = lambda pid: f"{_REACTOME_BASE}/data/pathways/low/entity/{pid}",
            parse_fn = lambda j: {p.get("stId", "").strip() for p in j if p.get("stId", "").strip()},
            n_workers = w,
            label    = "Reactome pathway",
        )

    def _reactome_proteins_by_pathways(self, pathway_ids: list[str]) -> dict[str, set[str]]:
        """Look up Reactome proteins for all R-* pathway IDs in parallel."""
        reactome_ids = [pid for pid in pathway_ids if self._is_reactome_id(pid)]
        result: dict[str, set[str]] = {pid: set() for pid in pathway_ids}
        if reactome_ids:
            w = max(1, min(len(reactome_ids), self._MAX_WORKERS))
            result.update(self._parallel_get(
                reactome_ids,
                url_fn   = lambda pid: f"{_REACTOME_BASE}/data/participants/{pid}/referenceEntities",
                parse_fn = lambda j: {
                    e.get("identifier", "").strip()
                    for e in j
                    if e.get("databaseName") == "UniProt" and e.get("identifier", "").strip()
                },
                n_workers = w,
                label    = "Reactome protein",
            ))
        return result

    # ── ChEMBL (chunked __in filter + pagination; parallel for reverse) ───────

    def _chembl_pathways_by_proteins(self, protein_ids: list[str]) -> dict[str, set[str]]:
        """Fetch ChEMBL pathway xrefs for all proteins.

        Chunks into groups of _CHEMBL_CHUNK to stay under URL length limits,
        then paginates each chunk. Only keeps xrefs from known pathway databases.
        """
        result: dict[str, set[str]] = {pid: set() for pid in protein_ids}
        # ChEMBL's __in filter has a hard URL length limit, so we batch the
        # accessions in chunks of _CHEMBL_CHUNK (default 50).
        for i in range(0, len(protein_ids), self._CHEMBL_CHUNK):
            chunk   = protein_ids[i : i + self._CHEMBL_CHUNK]
            targets = self._chembl_paginate(
                f"{_CHEMBL_BASE}/target"
                f"?target_components__accession__in={','.join(chunk)}&format=json&limit=100"
            )
            for target in targets:
                # Each target record contains a list of protein components.
                # We only want components whose accession is in our request set.
                for component in target.get("target_components", []):
                    acc = component.get("accession", "").strip()
                    if acc not in result:
                        continue
                    # Each component can have cross-references to external pathway DBs.
                    for xref in component.get("target_component_xrefs", []):
                        if xref.get("xref_src_db") in self._CHEMBL_PATHWAY_DBS:
                            xref_id = xref.get("xref_id", "").strip()
                            if xref_id:
                                result[acc].add(xref_id)
        return result

    def _chembl_proteins_by_pathways(self, pathway_ids: list[str]) -> dict[str, set[str]]:
        """Fetch ChEMBL proteins for KEGG/WikiPathways/PharmGKB IDs in parallel.

        ChEMBL returns 400 for SMP* and R-HSA-* IDs, so those are skipped.
        """
        chembl_ids = [
            pid for pid in pathway_ids
            if not self._is_smpdb_id(pid) and not self._is_reactome_id(pid)
        ]
        result: dict[str, set[str]] = {pid: set() for pid in pathway_ids}
        if chembl_ids:
            w = max(1, min(len(chembl_ids), self._MAX_WORKERS))
            result.update(self._parallel_get(
                chembl_ids,
                url_fn   = lambda pid: f"{_CHEMBL_BASE}/target_component?target_component_xrefs__xref_id={pid}&format=json&limit=200",
                parse_fn = lambda j: {c.get("accession", "").strip() for c in j if c.get("accession", "").strip()},
                n_workers = w,
                label    = "ChEMBL protein",
            ))
        return result

    # ── Public interface ──────────────────────────────────────────────────────

    def get_pathways_by_protein(self, protein_ids: str | list[str]) -> dict[str, list[str]]:
        """Return pathway IDs for each protein, merged across SMPDB, Reactome, ChEMBL,
        and optionally DrugBank.

        Args:
            protein_ids: A UniProt accession or list of accessions.

        Returns:
            Dict of {accession: sorted_list_of_pathway_ids}. Pathway IDs span
            multiple namespaces: SMP* (SMPDB), R-HSA-* (Reactome), hsa* (KEGG),
            WP* (WikiPathways), etc.

        Example:
            mapper = PathwayMapper()
            result = mapper.get_pathways_by_protein(["P00734", "P00533"])
        """
        ids = self._normalize(protein_ids)

        db       = self._drugbank_pathways_by_proteins(set(ids))
        smpdb    = self._smpdb_pathways_by_proteins(ids)
        reactome = self._reactome_pathways_by_proteins(ids)
        chembl   = self._chembl_pathways_by_proteins(ids)

        return {
            pid: sorted(db[pid] | smpdb[pid] | reactome[pid] | chembl[pid])
            for pid in ids
        }

    def get_proteins_by_pathway(self, pathway_ids: str | list[str]) -> dict[str, list[str]]:
        """Return UniProt accessions for all proteins in each pathway, merged across sources.

        Handles mixed pathway namespaces in a single call — pass SMPDB, Reactome, KEGG,
        and WikiPathways IDs together and get back a single merged result per pathway.

        Args:
            pathway_ids: A pathway ID or list of pathway IDs.

        Returns:
            Dict of {pathway_id: sorted_list_of_uniprot_accessions}.

        Example:
            mapper = PathwayMapper()
            result = mapper.get_proteins_by_pathway(["SMP0000278", "R-HSA-159740"])
        """
        ids = self._normalize(pathway_ids)

        db       = self._drugbank_proteins_by_pathways(set(ids))
        smpdb    = self._smpdb_proteins_by_pathways(ids)
        reactome = self._reactome_proteins_by_pathways(ids)
        chembl   = self._chembl_proteins_by_pathways(ids)

        return {
            pid: sorted(db[pid] | smpdb[pid] | reactome[pid] | chembl[pid])
            for pid in ids
        }

# ── Quick demo ────────────────────────────────────────────────────────────────

if __name__ == "__main__":
    import os
    SEP = "-" * 60

    # --- Four core framework operations ---

    print(SEP)
    print("Prot(s)  —  get_proteins_for_pathway('R-HSA-159740')")
    proteins = get_proteins_for_pathway("R-HSA-159740")
    print(f"  {len(proteins)} proteins: {proteins[:5]}{'...' if len(proteins) > 5 else ''}")

    print(SEP)
    print("Π(u)  —  get_pathways_for_protein('P00734')")
    pathways = get_pathways_for_protein("P00734")
    smpdb_pw = [p for p in pathways if p.startswith("SMP")]
    react_pw = [p for p in pathways if p.startswith("R-")]
    print(f"  {len(pathways)} pathways total  |  SMPDB={len(smpdb_pw)}  Reactome={len(react_pw)}")
    print(f"  First 5: {pathways[:5]}")

    print(SEP)
    print("Φ_native(d)  —  phi_native(['R-HSA-159740', 'R-HSA-140837'])")
    native = phi_native(["R-HSA-159740", "R-HSA-140837"])
    print(f"  {len(native)} proteins in union: {sorted(native)[:5]}{'...' if len(native) > 5 else ''}")

    print(SEP)
    print("Φ_infer(d; l)  —  phi_infer(['P00734', 'P00740'])  (two-hop expansion)")
    inferred = phi_infer(["P00734", "P00740"])
    print(f"  Pathways : {len(inferred['pathways'])} — first 5: {inferred['pathways'][:5]}")
    print(f"  Neighbors: {len(inferred['neighbors'])} proteins")

    print(SEP)
    print("phi_infer_chembl('CHEMBL1536', min_pchembl=6)  — warfarin via ChEMBL (MoA + bioactivity)")
    chembl_result = phi_infer_chembl("CHEMBL1536", min_pchembl=6)
    print(f"  ChEMBL proteins : {len(chembl_result['proteins'])} — {chembl_result['proteins'][:5]}")
    print(f"  Pathways        : {len(chembl_result['pathways'])} — first 5: {chembl_result['pathways'][:5]}")
    print(f"  Neighbors       : {len(chembl_result['neighbors'])} proteins")

    print(SEP)
    print("get_sequence('P00734')  — first 30 residues")
    seq = get_sequence("P00734")
    if seq:
        print(f"  {seq[:30]}...")

    print(SEP)
    print("get_go_terms('P00734')")
    go = get_go_terms("P00734")
    print(f"  {go[:5]}{'...' if len(go) > 5 else ''}")

    print(SEP)
    print("get_pfams('P00734')")
    pfam = get_pfams("P00734")
    print(f"  {pfam}")

    # --- PathwayMapper class demo (with optional DrugBank XML) ---
    xml_file = r"C:\Users\ashto\OneDrive\DDI codes\DrugBank Parsing Code\drugbank_all_full_database_V5.1.14.xml\drugbank_full_database_V5.1.14.xml"

    if os.path.exists(xml_file):
        print(SEP)
        print("PathwayMapper — get_proteins_by_pathway(['SMP0000278', 'R-HSA-159740'])  (with DrugBank)")
        mapper  = PathwayMapper(xml_file)
    else:
        print(SEP)
        print("PathwayMapper — get_proteins_by_pathway(['SMP0000278', 'R-HSA-159740'])  (no DrugBank XML)")
        mapper  = PathwayMapper()

    result = mapper.get_proteins_by_pathway(["SMP0000278", "R-HSA-159740"])
    for pid, prots in result.items():
        print(f"  {pid} ({len(prots)} proteins): {prots[:5]}{'...' if len(prots) > 5 else ''}")
