"""
Pulls pathway-to-protein and protein-to-pathway mappings out of the DrugBank XML.

The two functions here are the building blocks for the pathway profiles used
in DDI prediction:

  get_proteins_by_pathway  ->  given an SMPDB pathway ID, which proteins are in it?
  get_pathways_by_protein  ->  given a UniProt protein ID, which pathways does it appear in?

These map to Prot(s) and Pi(u) in the biological-profile framework.
The native pathway profile for a drug pools Prot(s) across all of its pathways.
The inferred profile does a two-hop walk: protein -> pathways -> all proteins in those pathways.
"""

import xml.etree.ElementTree as ET
import zipfile
import gzip
import time

# All DrugBank XML tags are namespaced under this URI
NS = "{http://www.drugbank.ca}"


def _open_xml(xml_path): # This might be able to become an importable utility function if we need to read other DrugBank XML files in the future
    """Returns a readable file handle for a DrugBank XML file (.xml, .zip, or .gz)."""
    if xml_path.endswith(".zip"):
        zf = zipfile.ZipFile(xml_path, "r")
        xml_files = [f for f in zf.namelist() if f.endswith(".xml")]
        if not xml_files:
            raise FileNotFoundError("No XML file found inside the zip archive.")
        return zf.open(xml_files[0])
    elif xml_path.endswith(".gz"):
        return gzip.open(xml_path, "rb")
    else:
        return open(xml_path, "rb")


def get_proteins_by_pathway(xml_path, pathway_id):
    """
    Given an SMPDB pathway ID, return the UniProt IDs of all proteins in that pathway.

    This is Prot(s) from the framework. DrugBank lists the member proteins of each
    pathway under <pathway>/<enzymes>/<uniprot-id>. The same pathway can show up
    under multiple drugs in the XML, so we union proteins across all occurrences.

    Used to build:
      - Phi_native: pool Prot(s) for every pathway in a drug's own <pathways> list
      - Phi_infer:  the second hop after Pi(u) returns which pathways a protein is in

    Example:
        get_proteins_by_pathway("drugbank.xml", "SMP00112")
        # -> ["P00734", "P00748", "P03952"]
    """
    proteins = set()

    with _open_xml(xml_path) as fh:
        # Iterate at the <drug> level so the full <pathways> subtree is still in
        # memory when we read it. Clearing at a finer grain (e.g. every non-pathway
        # element) destroys child nodes before the parent pathway end-event fires,
        # which is why findall() would return nothing.
        for _event, elem in ET.iterparse(fh, events=["end"]):
            if elem.tag != f"{NS}drug":
                continue

            for pathway in elem.findall(f"{NS}pathways/{NS}pathway"):
                smpdb_id = pathway.findtext(f"{NS}smpdb-id")
                if smpdb_id and smpdb_id.strip() == pathway_id:
                    for uid in pathway.findall(f"{NS}enzymes/{NS}uniprot-id"):
                        if uid.text and uid.text.strip():
                            proteins.add(uid.text.strip())

            elem.clear()  # free the whole drug subtree once we're done with it

    return sorted(proteins)


def get_pathways_by_protein(xml_path, protein_id):
    """
    Given a UniProt protein ID, return the SMPDB pathway IDs it appears in.

    This is Pi(u) from the framework. It is the first hop in the Phi_infer
    two-hop expansion:

        for each protein u the drug acts on:
            find Pi(u)  <-- this function
            then for each pathway s in Pi(u), call get_proteins_by_pathway(s)
            pool all resulting proteins

    Example:
        get_pathways_by_protein("drugbank.xml", "P00734")
        # -> ["SMP00112", "SMP00765"]
    """
    pathways = set()

    with _open_xml(xml_path) as fh:
        # Iterate at the <drug> level so the full <pathways> subtree is still in
        # memory when we read it (same fix as get_proteins_by_pathway).
        for _event, elem in ET.iterparse(fh, events=["end"]):
            if elem.tag != f"{NS}drug":
                continue

            for pathway in elem.findall(f"{NS}pathways/{NS}pathway"):
                smpdb_id = pathway.findtext(f"{NS}smpdb-id")
                if smpdb_id:
                    for uid in pathway.findall(f"{NS}enzymes/{NS}uniprot-id"):
                        if uid.text and uid.text.strip() == protein_id:
                            pathways.add(smpdb_id.strip())
                            break  # found it in this pathway, move on

            elem.clear()  # free the whole drug subtree once we're done with it

    return sorted(pathways)


def _print_first_pathway(xml_path):
    """Debug helper: print the raw XML of the first <pathway> found in the file.
    Run this to confirm the exact tag structure before querying."""
    with _open_xml(xml_path) as fh:
        for _event, elem in ET.iterparse(fh, events=["end"]):
            if elem.tag != f"{NS}drug":
                continue
            pathway = elem.find(f"{NS}pathways/{NS}pathway")
            if pathway is not None:
                print(ET.tostring(pathway, encoding="unicode"))
                break
            elem.clear()


if __name__ == "__main__":
    xml_file = r"C:\Users\ashto\OneDrive\DDI codes\DrugBank Parsing Code\drugbank_all_full_database_V5.1.14.xml\drugbank_full_database_V5.1.14.xml"

    # Test get_proteins_by_pathway: Lepirudin Action Pathway (confirmed in XML)
    print("Proteins in SMP0000278 (Lepirudin Action Pathway)")
    proteins = get_proteins_by_pathway(xml_file, "SMP0000278")
    print(proteins)

    # Test get_pathways_by_protein: P00734 (thrombin, confirmed in XML above)
    print("\nPathways for P00734 (thrombin)")
    pathways = get_pathways_by_protein(xml_file, "P00734")
    print(pathways)