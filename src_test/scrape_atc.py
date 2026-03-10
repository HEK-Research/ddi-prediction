"""
Scrapes the Anatomical-Therapeutic-Chemical (ATC) classes from the WHO website (https://https://atcddd.fhi.no/atc_ddd_index/).
It reads ATC classes and their information, and writes to one flat CSV file. 
Python conversion of atcd.R by Fabricio Kury https://github.com/fabkury/atcd/tree/master
"""
import os
import pickle
import requests
import pandas as pd
from bs4 import BeautifulSoup
from datetime import date

# ------- Globals -------
OUT_DIR = os.path.join(os.path.dirname(__file__),"..","data", "atc_code_output")
RDS_DIR = os.path.join(OUT_DIR, "rds_cache")
os.makedirs(OUT_DIR, exist_ok=True)
os.makedirs(RDS_DIR, exist_ok=True)

BASE_URL = "https://atcddd.fhi.no/atc_ddd_index"

ATC_ROOTS = ['A', 'B', 'C', 'D', 'G', 'H', 'J', 'L', 'M', 'N', 'P', 'R', 'S', 'V']

# ------- Caching Utilities -------
def _cache_path(name: str) -> str:
    return os.path.join(OUT_DIR, f"{name}.pkl")

def wrap_cache(name:str, build_fn):
    """Load from cache if available, otherwise call build_fn(), cache, and return the result."""
    path = _cache_path(name)
    if os.path.exists(path):
        print(f"Reading '{name}' from cache '{path}'...")
        with open(path, 'rb') as f:
            return pickle.load(f)
    print(f"Building '{name}'...")
    val = build_fn()
    print(f"'{name}' completed. Saving to '{path}'...")
    with open(path, 'wb') as f:
        pickle.dump(val, f)
    return val

def get_cache(name:str):
    """Load from cache or raise an error."""
    path = _cache_path(name)
    if not os.path.exists(path):
        raise FileNotFoundError(f"Cache '{name}' not found at '{path}'.")
    print(f"Reading '{name}' from cache '{path}'...")
    with open(path, 'rb') as f:
        return pickle.load(f)

# -------- Scraper --------
def scrape_who_atc(root_atc_code):
    """
    Recursively scrapes ATC codes for root_atc_code and its children,
    stopping at level 4. Names are intentionally not collected.
    Returns a flat list of ATC code strings (levels 1-4).
    """
    code_length = len(root_atc_code)

    if code_length >= 5:
        return []

    url = f"{BASE_URL}/?code={root_atc_code}&showdescription=no"
    print(f"Scraping {url}")

    response = requests.get(url, timeout=30)
    response.raise_for_status()
    soup = BeautifulSoup(response.text, "html.parser")

    codes = []

    if code_length == 1:
        codes.append(root_atc_code)

    # ATC code lengths by level: L1=1, L2=3, L3=4, L4=5, L5=7
    # Expected child length given current code length
    next_length = {1: 3, 3: 4, 4: 5}
    expected_child_length = next_length.get(code_length)
    if expected_child_length is None:
        return codes  # no children expected

    all_links = soup.find_all("a", href=lambda h: h and "code=" in h and "showdescription" in h)

    for link in all_links:
        href = link.get("href", "")
        child_code = href.split("code=")[1].split("&")[0].strip()
        child_length = len(child_code)

        # Only collect direct children at the expected next length
        if child_length != expected_child_length:
            continue

        codes.append(child_code)

        # Recurse if child is still at or below level 4 (length <= 5)
        if child_length <= 4:
            codes.extend(scrape_who_atc(child_code))

    return codes

# -------- Main --------
if __name__ == "__main__":
    for root in ATC_ROOTS:
        cache_name = f"who_atc_{root}"
        wrap_cache(cache_name, lambda r=root: scrape_who_atc(r))

    all_codes = []
    for root in ATC_ROOTS:
        all_codes.extend(get_cache(f"who_atc_{root}"))

    # Map code length to ATC level number and encoding weight
    length_to_level  = {1: 1, 3: 2, 4: 3, 5: 4}
    length_to_weight = {3: 1.0, 4: 2.0, 5: 3.0}  # L2=1, L3=2, L4=3

    # Discard Level 1 codes (length 1) — too coarse for encoding
    encoding_codes = [c for c in all_codes if len(c) > 1]

    # Sort: L2 first, then L3, then L4; alphabetically within each level
    encoding_codes.sort(key=lambda c: (len(c), c))

    who_atc = pd.DataFrame({
        "atc_code": encoding_codes,
        "level":    [length_to_level[len(c)]  for c in encoding_codes],
        "weight":   [length_to_weight[len(c)] for c in encoding_codes],
    })

    # The position (0-based row index) is the bit position in the encoding vector
    who_atc.insert(0, "index", range(len(who_atc)))

    out_file = os.path.join(OUT_DIR, f"WHO_ATC_codes_{date.today().isoformat()}.csv")
    print(f"Writing results to {out_file}")
    if os.path.exists(out_file):
        print("Warning: file already exists. Will be overwritten.")
    who_atc.to_csv(out_file, index=False)

    print(f"Script execution completed. Total encoding codes: {len(who_atc)}")
    print(who_atc.groupby("level").size().rename("count"))
        
