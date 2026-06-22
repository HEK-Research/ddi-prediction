"""
smpdb_protein_pathway_batch.py
================================================================================
 
WHAT THIS SCRIPT DOES
---------------------
Builds batch-level pathway-protein networks using:
 
1. SMPDB protein dump (canonical)        — local ZIP of pathway CSVs
2. DrugBank XML (optional enrichment)    — local ZIP/GZIP/XML file
3. Reactome API (live external lookups)  — REST calls to reactome.org
4. ChEMBL drug-target mappings           — REST calls to EBI ChEMBL API
5. UniProt feature retrieval             — REST calls to UniProt REST API
 
Supported query patterns:
- Pathway  → proteins      (direct lookup)
- Protein  → pathways      (reverse lookup)
- Two-hop pathway expansions (protein → pathway → neighbors)
- ChEMBL drug → protein neighborhood expansion
- UniProt annotations: sequence, GO terms, Pfam domains
 
PERFORMANCE NOTES
-----------------
- All external API calls use a shared requests.Session for TCP connection reuse.
- Reactome and UniProt calls are parallelised with ThreadPoolExecutor.
- SMPDB CSV loading skips iterrows() in favour of vectorised groupby.
- DrugBank XML is stream-parsed (iterparse) to avoid loading the full tree.
- ChEMBL target → UniProt resolution is batched in chunks of 50.
"""
 
import gzip
import xml.etree.ElementTree as ET
import zipfile
import os
from concurrent.futures import ThreadPoolExecutor, as_completed  # parallel HTTP requests
 
import pandas as pd
import requests
from requests.adapters import HTTPAdapter          # connection pool config
from urllib3.util.retry import Retry               # automatic retry on transient failures
 
# ─────────────────────────────────────────────────────────────
# Shared HTTP Session (reuse TCP connections across all calls)
# ─────────────────────────────────────────────────────────────
 
def _build_session(
    total_retries: int = 3,
    backoff_factor: float = 0.4,
    pool_connections: int = 10,
    pool_maxsize: int = 20
) -> requests.Session:
    """
    Create a requests.Session configured with:
      - Automatic retry on HTTP 429 / 500 / 502 / 503 / 504 errors.
      - Exponential back-off between retries (backoff_factor).
      - A connection pool large enough for parallel Reactome/UniProt calls.
 
    Using a Session instead of bare requests.get() reuses the underlying
    TCP connections (HTTP keep-alive), which meaningfully cuts latency when
    you are making dozens of calls to the same host.
    """
    session = requests.Session()
 
    # Retry policy: retry on the listed HTTP status codes, with exponential wait.
    # backoff_factor=0.4 → waits 0s, 0.4s, 0.8s, 1.6s … between attempts.
    retry = Retry(
        total=total_retries,
        backoff_factor=backoff_factor,
        status_forcelist=[429, 500, 502, 503, 504],
        allowed_methods=["GET"],       # only retry idempotent GETs
        raise_on_status=False          # we check status ourselves
    )
 
    # Mount the adapter on both http:// and https:// prefixes.
    adapter = HTTPAdapter(
        max_retries=retry,
        pool_connections=pool_connections,   # number of distinct hosts kept alive
        pool_maxsize=pool_maxsize            # max connections per host
    )
    session.mount("https://", adapter)
    session.mount("http://", adapter)
 
    return session
 
 
# Module-level shared session — all API helpers below import this.
_SESSION = _build_session()
 
# XML namespace prefix used throughout the DrugBank XML schema.
NS = "{http://www.drugbank.ca}"
 
# Base URLs for each external API — centralised so they are easy to update.
_REACTOME_BASE = "https://reactome.org/ContentService"   # Reactome REST API
_CHEMBL_BASE   = "https://www.ebi.ac.uk/chembl/api/data" # ChEMBL REST API
_UNIPROT_BASE  = "https://rest.uniprot.org"              # UniProt REST API
 
 
# ─────────────────────────────────────────────────────────────
# UniProt Batch Feature Extractor
# ─────────────────────────────────────────────────────────────
 
class UniprotConverter:
    """
    Retrieve sequence, GO terms, and Pfam domain annotations from UniProt
    for a list of accession IDs.
 
    A single `_fetch_batch()` call downloads all fields in one HTTP round-trip
    per 100 IDs, and the result is cached on the instance so that
    `uniprot_to_sequence()`, `uniprot_to_GO_terms()`, and `uniprot_to_pfams()`
    can all be called without repeating network I/O.
    """
 
    def __init__(self, uniprot_ids):
        # Accept either a single string or an iterable of accession strings.
        # dict.fromkeys preserves insertion order while deduplicating, which is
        # faster than list(set(...)) because it avoids hashing overhead on large lists.
        self.uniprot_ids = (
            [uniprot_ids]
            if isinstance(uniprot_ids, str)
            else list(dict.fromkeys(uniprot_ids))
        )
 
        # Lazy-load cache: None means "not yet fetched".
        # After the first fetch, this holds {accession: entry_dict}.
        self._cached_payload = None
 
    def _fetch_batch(self) -> dict:
        """
        Download UniProt entries for all accessions in self.uniprot_ids.
 
        Strategy:
        - Requests are sent in chunks of 100 (UniProt recommends ≤200 per query).
        - We request only the fields we actually use: accession, sequence, go,
          and Pfam cross-references. Fetching fewer fields speeds up the response.
        - Results are merged into a single {accession: raw_entry} dict and cached.
 
        Returns:
            dict mapping UniProt accession → raw JSON entry dict from UniProt API.
        """
 
        # Return the previously fetched result if available (avoid re-fetching).
        if self._cached_payload is not None:
            return self._cached_payload
 
        # Nothing to fetch — return an empty mapping immediately.
        if not self.uniprot_ids:
            return {}
 
        chunk_size = 100          # max IDs per UniProt query (safe limit)
        combined_results = {}     # accumulates {accession: entry} across all chunks
 
        for i in range(0, len(self.uniprot_ids), chunk_size):
            # Slice out the current chunk of up to 100 accession IDs.
            chunk = self.uniprot_ids[i:i + chunk_size]
 
            # Build a UniProt search query: "(accession:P00734) OR (accession:P00740) …"
            # Each accession is wrapped in parentheses as required by UniProt query syntax.
            query_string = " OR ".join(
                [f"(accession:{acc})" for acc in chunk]
            )
 
            # Request only the fields we need; omitting unnecessary fields reduces
            # payload size and speeds up both the server response and JSON parsing.
            params = {
                "query":  query_string,
                "fields": "accession,sequence,go,xref_pfam",
                "format": "json",
                "size":   len(chunk)   # cap the result count to exactly our chunk size
            }
 
            # Fire the GET request using the shared session (connection reuse).
            resp = _SESSION.get(
                f"{_UNIPROT_BASE}/uniprotkb/search",
                params=params,
                timeout=30   # seconds before giving up — prevents indefinite hangs
            )
 
            # Skip chunks that return an error rather than crashing the whole run.
            if resp.status_code != 200:
                print(f"[UniProt] HTTP {resp.status_code} for chunk {i}–{i+chunk_size}")
                continue
 
            # Parse the JSON body; "results" is a list of UniProt entry objects.
            results = resp.json().get("results", [])
 
            # Index each entry by its primary accession for O(1) lookup later.
            for entry in results:
                acc = entry.get("primaryAccession")
                if acc:
                    combined_results[acc] = entry
 
        # Store the merged result so subsequent method calls don't re-fetch.
        self._cached_payload = combined_results
        return self._cached_payload
 
    def uniprot_to_sequence(self) -> dict:
        """
        Return {accession: amino_acid_sequence_string} for all loaded IDs.
 
        The sequence lives at entry["sequence"]["value"] in the UniProt JSON schema.
        Missing sequences are returned as empty strings rather than raising KeyError.
        """
        batch = self._fetch_batch()   # uses cache if already fetched
 
        return {
            acc: data.get("sequence", {}).get("value", "")
            for acc, data in batch.items()
            # .get("sequence", {}) avoids KeyError if the field is absent
            # .get("value", "")    returns "" rather than None for missing sequences
        }
 
    def uniprot_to_GO_terms(self) -> dict:
        """
        Return {accession: [GO_id, …]} for all loaded IDs.
 
        GO terms are stored as cross-references in the UniProt JSON under
        uniProtKBCrossReferences, each with a "database" field equal to "GO".
        We filter those and collect just the "id" field (e.g. "GO:0006508").
        """
        batch = self._fetch_batch()
        go_map = {}
 
        for acc, data in batch.items():
            # uniProtKBCrossReferences is a list of {database, id, properties} dicts.
            xrefs = data.get("uniProtKBCrossReferences", [])
 
            # Keep only cross-references whose database is "GO" and extract the id.
            go_map[acc] = [
                x["id"]
                for x in xrefs
                if x.get("database") == "GO"
            ]
 
        return go_map
 
    def uniprot_to_pfams(self) -> dict:
        """
        Return {accession: [Pfam_id, …]} for all loaded IDs.
 
        Pfam entries are stored the same way as GO terms but with database == "Pfam".
        Pfam IDs look like "PF00051" (serine protease inhibitor domain, etc.).
        """
        batch = self._fetch_batch()
        pfam_map = {}
 
        for acc, data in batch.items():
            xrefs = data.get("uniProtKBCrossReferences", [])
 
            # Collect only Pfam cross-references by filtering on database name.
            pfam_map[acc] = [
                x["id"]
                for x in xrefs
                if x.get("database") == "Pfam"
            ]
 
        return pfam_map
 
 
# ─────────────────────────────────────────────────────────────
# ChEMBL Target Resolver
# ─────────────────────────────────────────────────────────────
 
def get_targets_from_chembl_batch(
    molecule_chembl_ids: list,
    min_pchembl: float = None
) -> dict:
    """
    Resolve a list of ChEMBL molecule IDs → UniProt protein accessions.
 
    Two complementary ChEMBL data sources are queried:
    1. /mechanism  — curated mechanism-of-action entries (high confidence).
    2. /activity   — measured bioactivity values filtered to SINGLE_PROTEIN targets
                     with a non-null pChEMBL value (optional minimum threshold).
 
    Both sets of target ChEMBL IDs are then resolved to UniProt accessions
    via /target endpoint, batched in chunks of 50.
 
    Args:
        molecule_chembl_ids: list of ChEMBL molecule IDs, e.g. ["CHEMBL1536"].
        min_pchembl: optional minimum pChEMBL value (log-scale affinity).
                     pChEMBL ≥ 6 corresponds roughly to IC50 ≤ 1 µM.
 
    Returns:
        dict mapping each molecule_id → sorted list of UniProt accession strings.
    """
 
    # Deduplicate the molecule list while preserving order.
    molecule_ids = list(dict.fromkeys(molecule_chembl_ids))
 
    # Pre-populate the output dict with empty sets so every input key is present
    # even if no targets are found for it.
    molecule_to_targets = {m_id: set() for m_id in molecule_ids}
 
    # Intermediate mapping: ChEMBL target ID → set of molecule IDs that hit it.
    # Used to connect activity records back to their source molecule.
    target_to_mol = {}
 
    def _paginate(url: str) -> list:
        """
        Walk through all pages of a ChEMBL paginated endpoint and return a
        flat list of all records.
 
        ChEMBL uses a page_meta.next field containing the relative URL of the
        next page, or null when the last page is reached.
        The key for the records varies by endpoint, so we detect it from the URL.
        """
        accumulated = []    # all records collected across pages
        current_url = url   # start with the first page URL
 
        while current_url:
            r = _SESSION.get(current_url, timeout=30)  # shared session reuses connection
 
            if r.status_code != 200:
                break   # stop pagination on error rather than raising
 
            data = r.json()
 
            # The records key name differs by ChEMBL endpoint:
            # /mechanism → "mechanisms", /activity → "activities", /target → "targets"
            if "mechanism" in url:
                accumulated.extend(data.get("mechanisms", []))
            elif "activity" in url:
                accumulated.extend(data.get("activities", []))
            else:
                accumulated.extend(data.get("targets", []))
 
            # Fetch the relative URL of the next page (None when done).
            next_page = data.get("page_meta", {}).get("next")
 
            # Prepend the ChEMBL base URL to turn the relative path into an absolute URL.
            current_url = (
                f"{_CHEMBL_BASE}{next_page}"
                if next_page else None   # None exits the while loop
            )
 
        return accumulated
 
    # ── Step 1: Query mechanism-of-action data ──────────────────────────────
    # __in= accepts a comma-separated list, letting us batch all molecules in
    # one request series rather than one request per molecule.
    mechs = _paginate(
        f"{_CHEMBL_BASE}/mechanism?"
        f"molecule_chembl_id__in={','.join(molecule_ids)}"
        f"&format=json&limit=1000"   # 1000 per page to minimise round-trips
    )
 
    for m in mechs:
        t_id = m.get("target_chembl_id")    # ChEMBL target identifier
        m_id = m.get("molecule_chembl_id")  # ChEMBL molecule identifier
 
        if t_id and m_id:
            # Map each target to the set of molecules that modulate it.
            target_to_mol.setdefault(t_id, set()).add(m_id)
 
    # ── Step 2: Query bioactivity data ──────────────────────────────────────
    # Filter to SINGLE PROTEIN targets with a measured pChEMBL value so we
    # only get direct, quantitative protein-ligand interactions.
    activities = _paginate(
        f"{_CHEMBL_BASE}/activity?"
        f"molecule_chembl_id__in={','.join(molecule_ids)}"
        f"&target_type=SINGLE+PROTEIN"       # exclude complexes, cell lines, etc.
        f"&pchembl_value__isnull=false"       # require a measured affinity value
        f"&format=json&limit=1000"
    )
 
    for a in activities:
        pval = a.get("pchembl_value")   # pChEMBL = -log10(IC50/Ki/Kd in molar)
 
        # Apply the optional minimum affinity threshold.
        if min_pchembl is not None:
            try:
                # Skip this record if pChEMBL is missing or below the threshold.
                if pval is None or float(pval) < min_pchembl:
                    continue
            except (TypeError, ValueError):
                # pval was non-numeric — skip rather than crash.
                continue
 
        t_id = a.get("target_chembl_id")
        m_id = a.get("molecule_chembl_id")
 
        if t_id and m_id:
            target_to_mol.setdefault(t_id, set()).add(m_id)
 
    # Early exit: no targets found from either source.
    if not target_to_mol:
        return {m_id: [] for m_id in molecule_ids}
 
    # ── Step 3: Resolve ChEMBL target IDs → UniProt accessions ─────────────
    # We batch the resolution in chunks of 50 to stay within URL length limits.
    all_target_ids = list(target_to_mol.keys())
 
    for i in range(0, len(all_target_ids), 50):
        chunk = all_target_ids[i:i + 50]   # up to 50 ChEMBL target IDs per request
 
        r = _SESSION.get(
            f"{_CHEMBL_BASE}/target?"
            f"target_chembl_id__in={','.join(chunk)}"
            f"&format=json&limit=100",
            timeout=30
        )
 
        if r.status_code != 200:
            continue   # skip this chunk on error; partial results are still usable
 
        for target in r.json().get("targets", []):
            t_id = target.get("target_chembl_id")
 
            # target_components contains the protein subunits with UniProt accessions.
            for component in target.get("target_components", []):
                acc = component.get("accession")   # UniProt accession, e.g. "P00734"
 
                if acc and t_id in target_to_mol:
                    # Map the UniProt accession back to every molecule that hit this target.
                    for m_id in target_to_mol[t_id]:
                        molecule_to_targets[m_id].add(acc)
 
    # Convert sets to sorted lists for deterministic output and JSON-serializability.
    return {
        m_id: sorted(list(acc_set))
        for m_id, acc_set in molecule_to_targets.items()
    }
 
 
# ─────────────────────────────────────────────────────────────
# Pathway Mapper
# ─────────────────────────────────────────────────────────────
 
class PathwayMapper:
    """
    Central index for pathway ↔ protein relationships.
 
    Data is loaded from two optional local sources at construction time:
      - SMPDB ZIP: authoritative, fully offline, fast.
      - DrugBank XML: supplemental, stream-parsed to keep memory usage low.
 
    Live Reactome lookups are performed at query time to fill gaps.
 
    Internal data structures:
      pathway_to_proteins_map: {pathway_id: set(uniprot_id, …)}
      protein_to_pathways_map: {uniprot_id: set(pathway_id, …)}
 
    Both dicts are populated by every data source so that lookups in either
    direction are O(1).
    """
 
    def __init__(
        self,
        xml_path=None,              # path to DrugBank XML / ZIP / GZ file
        smpdb_protein_zip=None      # path to SMPDB proteins ZIP archive
    ):
        self.xml_path          = xml_path
        self.smpdb_protein_zip = smpdb_protein_zip
 
        # Both maps start empty; _load_smpdb_csv() and _parse_drugbank_xml()
        # populate them by merging data from their respective sources.
        self.pathway_to_proteins_map = {}   # {pathway_id: set(uniprot_id)}
        self.protein_to_pathways_map = {}   # {uniprot_id: set(pathway_id)}
 
        # ── Load SMPDB canonical data (if the file exists) ───────────────────
        # SMPDB is the primary source; DrugBank provides supplemental coverage.
        if smpdb_protein_zip and os.path.exists(smpdb_protein_zip):
            self._load_smpdb_csv()
 
        # ── Optionally enrich with DrugBank data ─────────────────────────────
        if xml_path and os.path.exists(xml_path):
            self._parse_drugbank_xml()
 
    def _open_xml(self):
        """
        Open the DrugBank XML regardless of the archive format.
 
        DrugBank distributes its database as:
          - A raw .xml file (development / small subset)
          - A .xml.gz GZIP archive
          - A .zip archive containing the XML
 
        Returns a file-like object suitable for xml.etree.ElementTree.iterparse().
        The caller is responsible for closing it (use as a context manager).
        """
 
        if not self.xml_path:
            raise ValueError("No DrugBank XML path was provided to PathwayMapper.")
 
        # ── Case 1: ZIP archive ───────────────────────────────────────────────
        if self.xml_path.endswith(".zip"):
            zf = zipfile.ZipFile(self.xml_path, "r")  # opened but not yet yielded
 
            # Find all .xml members inside the ZIP; take the first one.
            xml_candidates = [
                f for f in zf.namelist()
                if f.lower().endswith(".xml")
            ]
 
            if not xml_candidates:
                raise FileNotFoundError(
                    "No XML file found inside DrugBank ZIP archive."
                )
 
            selected_xml = xml_candidates[0]   # typically "full_database.xml"
            print(f"[DrugBank] ZIP detected → loading: {selected_xml}")
 
            # Return a file-like object for the inner XML without extracting it.
            return zf.open(selected_xml)
 
        # ── Case 2: GZIP archive ──────────────────────────────────────────────
        elif self.xml_path.endswith(".gz"):
            print("[DrugBank] GZIP detected")
            # gzip.open() returns a file-like object that decompresses on-the-fly.
            return gzip.open(self.xml_path, "rb")
 
        # ── Case 3: Raw uncompressed XML ──────────────────────────────────────
        elif self.xml_path.endswith(".xml"):
            print("[DrugBank] Raw XML detected")
            return open(self.xml_path, "rb")  # binary mode required by iterparse
 
        else:
            raise ValueError(
                f"Unsupported DrugBank file format: {self.xml_path}. "
                "Expected .xml, .xml.gz, or .zip"
            )
 
    def _load_smpdb_csv(self):
        """
        Parse all pathway-protein CSV files bundled inside the SMPDB ZIP.
 
        SMPDB distributes one CSV per pathway; each CSV has at minimum:
          - smpdb_id    : e.g. "SMP0000278"
          - uniprot_id  : e.g. "P00734"
 
        Performance: uses pandas vectorised groupby instead of iterrows(),
        which is typically 50–100× faster for large CSVs.
        """
 
        try:
            total_rows = 0   # running count of valid pathway-protein pairs loaded
 
            with zipfile.ZipFile(self.smpdb_protein_zip, "r") as zf:
 
                # Collect the names of all CSV members in the archive.
                csv_files = [
                    f for f in zf.namelist()
                    if f.endswith(".csv")
                ]
 
                print(f"[SMPDB] Found {len(csv_files)} pathway CSV files")
 
                for csv_name in csv_files:
 
                    with zf.open(csv_name) as f:
                        # Read the CSV; columns are normalised below.
                        df = pd.read_csv(f)
 
                    # Normalise column names: lowercase, strip whitespace, replace
                    # spaces with underscores so "SMPDB ID" → "smpdb_id".
                    df.columns = [
                        c.strip().lower().replace(" ", "_")
                        for c in df.columns
                    ]
 
                    # Skip files that don't have the two required columns.
                    if "smpdb_id" not in df.columns or "uniprot_id" not in df.columns:
                        print(f"[Warning] Skipping malformed file: {csv_name}")
                        continue
 
                    # Drop rows where either key column is missing or literally "nan".
                    df = df[["smpdb_id", "uniprot_id"]].dropna()
                    df = df[df["uniprot_id"].str.lower() != "nan"]
 
                    # Strip surrounding whitespace from both key columns.
                    df["smpdb_id"]   = df["smpdb_id"].str.strip()
                    df["uniprot_id"] = df["uniprot_id"].str.strip()
 
                    # ── Vectorised population of pathway_to_proteins_map ──────
                    # groupby smpdb_id, collect all uniprot_ids into a set per pathway.
                    for pathway_id, grp in df.groupby("smpdb_id"):
                        proteins = set(grp["uniprot_id"].tolist())
                        self.pathway_to_proteins_map.setdefault(
                            pathway_id, set()
                        ).update(proteins)
 
                    # ── Vectorised population of protein_to_pathways_map ──────
                    # Same idea reversed: groupby uniprot_id, collect pathway IDs.
                    for uniprot_id, grp in df.groupby("uniprot_id"):
                        pathways = set(grp["smpdb_id"].tolist())
                        self.protein_to_pathways_map.setdefault(
                            uniprot_id, set()
                        ).update(pathways)
 
                    total_rows += len(df)
 
            print(f"[SMPDB] Loaded {len(self.pathway_to_proteins_map)} pathways")
            print(f"[SMPDB] Loaded {len(self.protein_to_pathways_map)} proteins")
            print(f"[SMPDB] Loaded {total_rows} pathway-protein edges")
 
        except Exception as e:
            # Non-fatal: log and continue so DrugBank/Reactome can still work.
            print(f"[Warning] SMPDB load failed: {e}")
 
    def _parse_drugbank_xml(self):
        """
        Stream-parse DrugBank XML to extract supplemental pathway-protein edges.
 
        Uses ElementTree.iterparse() so that only one XML element is in memory
        at a time — essential because the full DrugBank XML exceeds 1 GB.
 
        Parsing strategy:
          - Accumulate SMPDB pathway IDs from <pathway><smpdb-id> elements.
          - Accumulate UniProt IDs from <polypeptide><external-identifiers> elements.
          - On </drug>, flush both sets into the maps and reset accumulators.
          - Call elem.clear() and root.clear() aggressively to free parsed memory.
        """
 
        try:
            with self._open_xml() as f:
 
                # iterparse yields (event, element) pairs without loading the full tree.
                context = ET.iterparse(f, events=("start", "end"))
                context = iter(context)
 
                # The first event is always the root element opening tag.
                event, root = next(context)
 
                # Accumulators for the drug currently being parsed.
                current_pathways    = set()   # SMPDB pathway IDs for this drug
                current_uniprot_ids = set()   # UniProt IDs for this drug's targets
                in_polypeptide      = False   # flag: are we inside a <polypeptide> block?
 
                for event, elem in context:
 
                    # ── Pathway ID extraction ─────────────────────────────────
                    # <pathway><smpdb-id>SMP0000001</smpdb-id></pathway>
                    if event == "end" and elem.tag == f"{NS}pathway":
                        smpdb_id = elem.find(f"{NS}smpdb-id")
 
                        if smpdb_id is not None and smpdb_id.text:
                            current_pathways.add(smpdb_id.text.strip())
 
                        # Free this element's memory immediately after reading.
                        elem.clear()
 
                    # ── Enter polypeptide block ───────────────────────────────
                    # We only want UniProt IDs from targets, which are under <polypeptide>.
                    elif event == "start" and elem.tag == f"{NS}polypeptide":
                        in_polypeptide = True
 
                    # ── UniProt ID extraction (inside polypeptide only) ───────
                    # <external-identifier><resource>UniProtKB</resource>
                    #   <identifier>P00734</identifier></external-identifier>
                    elif (
                        event == "end"
                        and elem.tag == f"{NS}external-identifier"
                        and in_polypeptide   # guard: only collect from target proteins
                    ):
                        resource   = elem.find(f"{NS}resource")
                        identifier = elem.find(f"{NS}identifier")
 
                        if (
                            resource   is not None
                            and identifier is not None
                            and resource.text == "UniProtKB"   # filter to UniProt only
                        ):
                            current_uniprot_ids.add(identifier.text.strip())
 
                        elem.clear()
 
                    # ── Exit polypeptide block ────────────────────────────────
                    elif event == "end" and elem.tag == f"{NS}polypeptide":
                        in_polypeptide = False
                        elem.clear()
 
                    # ── End of drug: flush accumulators into both maps ────────
                    elif event == "end" and elem.tag == f"{NS}drug":
 
                        # For each pathway found in this drug entry, add all the
                        # protein targets we collected from its polypeptide blocks.
                        for pid in current_pathways:
                            self.pathway_to_proteins_map.setdefault(
                                pid, set()
                            ).update(current_uniprot_ids)
 
                        # Symmetric reverse mapping: each protein → all its pathways.
                        for uid in current_uniprot_ids:
                            self.protein_to_pathways_map.setdefault(
                                uid, set()
                            ).update(current_pathways)
 
                        # Reset accumulators for the next drug record.
                        current_pathways    = set()
                        current_uniprot_ids = set()
 
                        # Clear the root element's children to free memory.
                        # Without this the entire tree accumulates in RAM.
                        root.clear()
 
        except Exception as e:
            # Non-fatal: log and continue; SMPDB data is still usable.
            print(f"[Warning] DrugBank parsing failed: {e}")
 
    def get_proteins_by_pathway_batch(self, target_pathways: list) -> dict:
        """
        Return the protein members of each requested pathway.
 
        For SMPDB and DrugBank pathways (SMP* prefix): served from the in-memory
        maps loaded at construction time — no network call needed.
 
        For Reactome pathways (R-HSA-* prefix): supplemented by a live call to
        the Reactome ContentService /participants endpoint.
 
        Reactome calls are issued in parallel using ThreadPoolExecutor to hide
        the per-request latency when querying many pathways at once.
 
        Args:
            target_pathways: list of pathway IDs (SMPDB or Reactome format).
 
        Returns:
            dict mapping each pathway_id → sorted list of UniProt accession strings.
        """
 
        # Seed results from the local in-memory map (O(1) per lookup).
        results = {
            pid: set(self.pathway_to_proteins_map.get(pid, set()))
            for pid in target_pathways
        }
 
        # Identify pathways that need a live Reactome lookup.
        reactome_ids = [pid for pid in target_pathways if "R-HSA" in pid]
 
        if reactome_ids:
            # ── Parallel Reactome fetch ───────────────────────────────────────
            def _fetch_reactome_proteins(rid: str):
                """Fetch UniProt IDs for one Reactome pathway; returns (rid, set)."""
                try:
                    resp = _SESSION.get(
                        f"{_REACTOME_BASE}/data/participants/"
                        f"{rid}/referenceEntities",
                        timeout=20
                    )
                    if resp.status_code == 200:
                        # referenceEntities contains proteins, metabolites, etc.
                        # We filter to UniProt-backed entries only.
                        return rid, {
                            entry["identifier"]
                            for entry in resp.json()
                            if entry.get("databaseName") == "UniProt"
                            and entry.get("identifier")
                        }
                except Exception:
                    pass
                return rid, set()   # return empty set on any failure
 
            # Run all Reactome fetches concurrently (I/O-bound, so threads help).
            # max_workers=10 matches the session's pool_connections setting.
            with ThreadPoolExecutor(max_workers=10) as executor:
                futures = {
                    executor.submit(_fetch_reactome_proteins, rid): rid
                    for rid in reactome_ids
                }
                for future in as_completed(futures):
                    rid, proteins = future.result()
                    results[rid].update(proteins)   # merge into the local result
 
        # Convert sets to sorted lists for deterministic output.
        return {k: sorted(list(v)) for k, v in results.items()}
 
    def get_pathways_by_protein_batch(self, target_proteins: list) -> dict:
        """
        Return all pathway memberships for each requested protein.
 
        Local data (SMPDB + DrugBank) is served from memory; Reactome pathways
        are fetched live in parallel.
 
        Args:
            target_proteins: list of UniProt accession strings.
 
        Returns:
            dict mapping each uniprot_id → sorted list of pathway ID strings.
        """
 
        # Seed from local maps.
        results = {
            p: set(self.protein_to_pathways_map.get(p, set()))
            for p in target_proteins
        }
 
        # ── Parallel Reactome reverse-lookup ─────────────────────────────────
        def _fetch_reactome_pathways(prot: str):
            """Fetch Reactome pathway IDs for one UniProt accession."""
            try:
                resp = _SESSION.get(
                    f"{_REACTOME_BASE}/data/mapping/UniProt/{prot}/pathways",
                    timeout=20
                )
                if resp.status_code == 200:
                    # Each element has a "stId" field like "R-HSA-159740".
                    return prot, {
                        path_node["stId"]
                        for path_node in resp.json()
                        if path_node.get("stId")
                    }
            except Exception:
                pass
            return prot, set()
 
        with ThreadPoolExecutor(max_workers=10) as executor:
            futures = {
                executor.submit(_fetch_reactome_pathways, prot): prot
                for prot in target_proteins
            }
            for future in as_completed(futures):
                prot, pathways = future.result()
                results[prot].update(pathways)
 
        return {k: sorted(list(v)) for k, v in results.items()}
 
 
# ─────────────────────────────────────────────────────────────
# Graph Expansion helpers
# ─────────────────────────────────────────────────────────────
 
def phi_native_batch(
    pathway_ids: list,
    xml_path=None,
    smpdb_protein_zip=None
) -> dict:
    """
    Direct (one-hop) pathway → protein expansion.
 
    Constructs a PathwayMapper from the provided data sources and returns a
    dict mapping each pathway ID to the set of UniProt accessions it contains.
 
    Args:
        pathway_ids:        list of SMPDB or Reactome pathway IDs.
        xml_path:           optional path to DrugBank XML / ZIP / GZ.
        smpdb_protein_zip:  optional path to SMPDB proteins ZIP archive.
 
    Returns:
        dict {pathway_id: set(uniprot_id, …)}
    """
 
    # Build the mapper (loads SMPDB and/or DrugBank from disk if paths given).
    mapper = PathwayMapper(
        xml_path=xml_path,
        smpdb_protein_zip=smpdb_protein_zip
    )
 
    # Query all pathways at once; convert list values to sets for downstream use.
    return {
        k: set(v)
        for k, v in mapper.get_proteins_by_pathway_batch(pathway_ids).items()
    }
 
 
def phi_infer_batch(
    protein_ids: list,
    xml_path=None,
    smpdb_protein_zip=None
) -> dict:
    """
    Two-hop protein neighborhood expansion.
 
    Hop 1: protein → all pathways it participates in.
    Hop 2: those pathways → all other proteins in them (neighbors).
 
    This is the core "guilt-by-association" step used in network pharmacology:
    a drug target protein is likely to functionally interact with other proteins
    that share its pathways.
 
    Args:
        protein_ids:        list of UniProt accession strings (seed proteins).
        xml_path:           optional path to DrugBank XML.
        smpdb_protein_zip:  optional path to SMPDB proteins ZIP.
 
    Returns:
        dict with keys:
          "pathways"  → sorted list of all discovered pathway IDs
          "neighbors" → sorted list of all co-pathway proteins (UniProt IDs)
    """
 
    mapper = PathwayMapper(
        xml_path=xml_path,
        smpdb_protein_zip=smpdb_protein_zip
    )
 
    # ── Hop 1: Seed proteins → their pathway memberships ─────────────────────
    pathway_map = mapper.get_pathways_by_protein_batch(protein_ids)
 
    # Flatten all discovered pathway IDs into a single deduplicated set.
    all_pathways = set()
    for pws in pathway_map.values():
        all_pathways.update(pws)
 
    # ── Hop 2: Discovered pathways → co-membership proteins ─────────────────
    protein_map = mapper.get_proteins_by_pathway_batch(list(all_pathways))
 
    # Flatten all protein neighbors from all pathways into one deduplicated set.
    all_neighbors = set()
    for prots in protein_map.values():
        all_neighbors.update(prots)
 
    return {
        "pathways":  sorted(list(all_pathways)),   # all pathways the seeds belong to
        "neighbors": sorted(list(all_neighbors))   # all proteins sharing those pathways
    }
 
 
def phi_infer_chembl_batch(
    molecule_chembl_ids: list,
    min_pchembl: float = None,
    xml_path=None,
    smpdb_protein_zip=None
) -> dict:
    """
    Full ChEMBL drug → pathway neighborhood expansion.
 
    Pipeline:
      1. Resolve each ChEMBL molecule → direct protein targets (via ChEMBL API).
      2. Map those targets → pathways (local SMPDB/DrugBank + live Reactome).
      3. Expand each pathway → all member proteins (the "neighborhood").
 
    This produces, for each drug:
      - proteins:   the direct molecular targets
      - pathways:   the pathways those targets participate in
      - neighbors:  all proteins co-present in those pathways
 
    Args:
        molecule_chembl_ids:  list of ChEMBL molecule IDs.
        min_pchembl:          optional minimum pChEMBL affinity threshold.
        xml_path:             optional DrugBank XML path.
        smpdb_protein_zip:    optional SMPDB ZIP path.
 
    Returns:
        dict mapping each molecule_id → {proteins, pathways, neighbors}.
    """
 
    # ── Step 1: Resolve molecule IDs → direct target UniProt accessions ──────
    # This makes ChEMBL mechanism + activity API calls (paginated).
    drug_to_targets = get_targets_from_chembl_batch(
        molecule_chembl_ids,
        min_pchembl=min_pchembl
    )
 
    # Collect the union of all target proteins across all drugs.
    # We query pathways for the whole set at once to avoid redundant API calls.
    all_proteins = set()
    for targets in drug_to_targets.values():
        all_proteins.update(targets)
 
    # ── Step 2: Global two-hop expansion (shared across all drugs) ───────────
    # Rather than running phi_infer_batch per drug (which would repeat the
    # Reactome calls), we run it once on the union of all targets.
    global_expansion = phi_infer_batch(
        list(all_proteins),
        xml_path=xml_path,
        smpdb_protein_zip=smpdb_protein_zip
    )
 
    # Reuse the same mapper instance for the per-drug breakdown queries below.
    mapper = PathwayMapper(
        xml_path=xml_path,
        smpdb_protein_zip=smpdb_protein_zip
    )
 
    # Build protein→pathways and pathway→proteins maps from the global expansion.
    protein_to_pathways = mapper.get_pathways_by_protein_batch(
        list(all_proteins)
    )
 
    pathway_to_proteins = mapper.get_proteins_by_pathway_batch(
        global_expansion["pathways"]   # only the pathways we already know about
    )
 
    # ── Step 3: Per-drug breakdown ────────────────────────────────────────────
    # Now split the global results back down to individual drug-level outputs.
    batch_output = {}
 
    for m_id in molecule_chembl_ids:
        drug_targets = drug_to_targets.get(m_id, [])   # direct protein targets
 
        # Collect pathways that at least one direct target belongs to.
        drug_pathways = set()
        for t in drug_targets:
            drug_pathways.update(protein_to_pathways.get(t, []))
 
        # Collect proteins that share any of those pathways (the neighborhood).
        drug_neighbors = set()
        for p in drug_pathways:
            drug_neighbors.update(pathway_to_proteins.get(p, []))
 
        # Store the three tiers of network context for this drug.
        batch_output[m_id] = {
            "proteins":  drug_targets,                    # direct ChEMBL targets
            "pathways":  sorted(list(drug_pathways)),     # implicated pathways
            "neighbors": sorted(list(drug_neighbors))     # co-pathway proteins
        }
 
    return batch_output
 
 
# ─────────────────────────────────────────────────────────────
# Demo / integration tests
# ─────────────────────────────────────────────────────────────
 
if __name__ == "__main__":
 
    # ── File paths (update these to your local copies) ───────────────────────
    smpdb_zip = r"C:\Users\ashto\ddi-prediction\data\smpdb_pathways_data_csv\smpdb_proteins.csv.zip"
    xml_file  = r"C:\Users\ashto\OneDrive\DDI codes\DrugBank Parsing Code\drugbank_all_full_database_V5.1.14.xml\drugbank_full_database_V5.1.14.zip"
 
    print("=" * 80)
    print("INITIALIZING PATHWAY MAPPER")
    print("=" * 80)
 
    # Create a shared mapper instance; this loads SMPDB CSVs and optionally
    # stream-parses DrugBank XML — can take a minute on large files.
    mapper = PathwayMapper(
        xml_path=xml_file,
        smpdb_protein_zip=smpdb_zip
    )
 
    print("\n")
 
    # ─────────────────────────────────────────────────────────
    # TEST 1: Direct pathway → protein lookup
    # ─────────────────────────────────────────────────────────
    print("=" * 80)
    print("TEST 1 — get_proteins_by_pathway_batch()")
    print("=" * 80)
 
    pathway_batch = mapper.get_proteins_by_pathway_batch([
        "SMP0000278",   # SMPDB: Blood Clotting Cascade
        "R-HSA-159740"  # Reactome: Hemostasis
    ])
 
    for pid, prots in pathway_batch.items():
        print(f"{pid} ({len(prots)} proteins)")
        print(prots[:10], "...\n")   # show first 10 to keep output readable
 
    # ─────────────────────────────────────────────────────────
    # TEST 2: Protein → pathway reverse lookup
    # ─────────────────────────────────────────────────────────
    print("=" * 80)
    print("TEST 2 — get_pathways_by_protein_batch()")
    print("=" * 80)
 
    protein_batch = mapper.get_pathways_by_protein_batch([
        "P00734",   # Prothrombin (F2)
        "P00740"    # Coagulation factor IX (F9)
    ])
 
    for prot, paths in protein_batch.items():
        # Split pathway IDs by source for a quick coverage summary.
        smp   = [p for p in paths if p.startswith("SMP")]
        react = [p for p in paths if p.startswith("R-HSA")]
 
        print(f"{prot}")
        print(f"  Total Pathways : {len(paths)}")
        print(f"  SMPDB          : {len(smp)}")
        print(f"  Reactome       : {len(react)}")
        print(f"  First 10       : {paths[:10]}\n")
 
    # ─────────────────────────────────────────────────────────
    # TEST 3: Native (one-hop) pathway expansion
    # ─────────────────────────────────────────────────────────
    print("=" * 80)
    print("TEST 3 — phi_native_batch()")
    print("=" * 80)
 
    native = phi_native_batch(
        ["SMP0000278", "R-HSA-159740"],
        xml_path=xml_file,
        smpdb_protein_zip=smpdb_zip
    )
 
    for pid, prots in native.items():
        print(f"{pid} ({len(prots)} proteins)")
        print(sorted(list(prots))[:10], "...\n")
 
    # ─────────────────────────────────────────────────────────
    # TEST 4: Two-hop protein expansion
    # ─────────────────────────────────────────────────────────
    print("=" * 80)
    print("TEST 4 — phi_infer_batch()")
    print("=" * 80)
 
    inferred = phi_infer_batch(
        ["P00734", "P00740"],
        xml_path=xml_file,
        smpdb_protein_zip=smpdb_zip
    )
 
    print(f"Discovered Pathways : {len(inferred['pathways'])}")
    print(f"First 10 pathways   : {inferred['pathways'][:10]}")
    print(f"Neighbor Proteins   : {len(inferred['neighbors'])}")
    print(f"First 20 neighbors  : {inferred['neighbors'][:20]}")
    print()
 
    # ─────────────────────────────────────────────────────────
    # TEST 5: ChEMBL drug target expansion
    # ─────────────────────────────────────────────────────────
    print("=" * 80)
    print("TEST 5 — phi_infer_chembl_batch()")
    print("=" * 80)
 
    chembl_expansion = phi_infer_chembl_batch(
        [
            "CHEMBL1536",   # Warfarin — vitamin K antagonist / coumarin anticoagulant
            "CHEMBL43"      # Heparin  — polysaccharide anticoagulant
        ],
        min_pchembl=6.0,        # only include targets with IC50 ≤ 1 µM
        xml_path=xml_file,
        smpdb_protein_zip=smpdb_zip
    )
 
    for drug, payload in chembl_expansion.items():
        print(f"{drug}")
        print(f"  Direct Targets : {len(payload['proteins'])}")
        print(f"  Pathways       : {len(payload['pathways'])}")
        print(f"  Neighbors      : {len(payload['neighbors'])}")
        print(f"  First Targets  : {payload['proteins'][:10]}")
        print()
 
    # ─────────────────────────────────────────────────────────
    # TEST 6: UniProt sequence retrieval
    # ─────────────────────────────────────────────────────────
    print("=" * 80)
    print("TEST 6 — uniprot_to_sequence()")
    print("=" * 80)
 
    converter = UniprotConverter([
        "P00734",   # Prothrombin
        "P00740"    # Factor IX
    ])
 
    seqs = converter.uniprot_to_sequence()
 
    for uid, seq in seqs.items():
        print(f"{uid}: length={len(seq)}")
        print(seq[:60] + "...\n")   # show first 60 amino acids only
 
    # ─────────────────────────────────────────────────────────
    # TEST 7: UniProt GO terms
    # ─────────────────────────────────────────────────────────
    print("=" * 80)
    print("TEST 7 — uniprot_to_GO_terms()")
    print("=" * 80)
 
    go_terms = converter.uniprot_to_GO_terms()
 
    for uid, gos in go_terms.items():
        print(f"{uid}: {len(gos)} GO terms")
        print(gos[:10], "\n")   # print first 10 GO IDs
 
    # ─────────────────────────────────────────────────────────
    # TEST 8: UniProt Pfam domain annotations
    # ─────────────────────────────────────────────────────────
    print("=" * 80)
    print("TEST 8 — uniprot_to_pfams()")
    print("=" * 80)
 
    pfams = converter.uniprot_to_pfams()
 
    for uid, pf in pfams.items():
        print(f"{uid}: {len(pf)} Pfam domains")
        print(pf, "\n")
 
    print("=" * 80)
    print("ALL TESTS COMPLETE")
    print("=" * 80)