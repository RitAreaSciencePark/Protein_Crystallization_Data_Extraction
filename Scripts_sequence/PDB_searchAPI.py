import requests
import json
import gemmi
import csv
import os
import re
from concurrent.futures import ThreadPoolExecutor, as_completed

# Setup cache directory for PDB files
CACHE_DIR = os.path.join(os.path.dirname(os.path.abspath(__file__)), ".pdb_cache")
os.makedirs(CACHE_DIR, exist_ok=True)

def get_cached_pubmed_id(pdb_id):
    """Retrieve PubMed ID with caching."""
    cache_file = os.path.join(CACHE_DIR, f"{pdb_id.lower()}_pubmed.json")
    if os.path.exists(cache_file):
        with open(cache_file, "r") as f:
            return json.load(f).get("pubmed_id", "NA")
    
    try:
        entry_data = requests.get(
            f"https://data.rcsb.org/rest/v1/core/entry/{pdb_id}",
            timeout=10
        ).json()
        pubmed = entry_data.get("rcsb_primary_citation", {}).get("pdbx_database_id_pub_med", "NA")
        
        # Cache the result
        with open(cache_file, "w") as f:
            json.dump({"pubmed_id": pubmed}, f)
        return pubmed
    except Exception as e:
        print(f"⚠ Warning fetching PubMed ID for {pdb_id}: {e}")
        return "NA"

def get_cached_mmcif(pdb_id):
    """Retrieve mmCIF file with caching."""
    cache_file = os.path.join(CACHE_DIR, f"{pdb_id.lower()}.cif")
    if os.path.exists(cache_file):
        with open(cache_file, "r") as f:
            return f.read()
    
    try:
        url = f"https://files.rcsb.org/view/{pdb_id}.cif"
        response = requests.get(url, timeout=10)
        response.raise_for_status()
        content = response.text
        
        # Cache the result
        with open(cache_file, "w") as f:
            f.write(content)
        return content
    except Exception as e:
        print(f"⚠ Warning fetching mmCIF for {pdb_id}: {e}")
        return None

def get_ph_from_mmcif_or_details(block):
    raw_ph = block.find_value("_exptl_crystal_grow.pH")
    try:
        if raw_ph not in (None, "", ".", "?", "NA"):
            return float(raw_ph)
    except ValueError:
        pass
    details = block.find_value("_exptl_crystal_grow.pdbx_details")
    if not details:
        return None
    match = re.search(r"\bpH\s*[:=]?\s*([0-9]+(?:\.[0-9]+)?)", details, re.IGNORECASE)
    if match:
        try:
            return float(match.group(1))
        except ValueError:
            return None
    return None

# ....get_method_from_mmcif_or_details...
def get_method_from_mmcif_or_details(block):
    method = block.find_value("_exptl_crystal_grow.method")
    if method not in (None, "", ".", "?", "NA"):
        return method.strip()
    details = block.find_value("_exptl_crystal_grow.pdbx_details")
    if not details:
        return None
    match = re.search(r"\b(method|technique|approach)\s*[:=]?\s*([a-zA-Z0-9 \-]+)", details, re.IGNORECASE)
    if match:
        return match.group(2).strip()
    return None

#.....get_temperature_from_mmcif_or_details.....
def get_temperature_from_mmcif_or_details(block):
   
    # ---- 1. Try mmCIF temperature fields first ----
    temp = block.find_value("_exptl_crystal_grow.temp") 

    if temp not in (None, "", ".", "?", "NA"):
        try:
            return float(temp)
        except ValueError:
            pass  # fall back to pdbx_details

    # ---- 2. Fallback: extract from pdbx_details ----
    details = block.find_value("_exptl_crystal_grow.pdbx_details")
    if not details:
        return None

    # Match ONLY Kelvin values like 277K or 277.00K
    match = re.search(
        r"\b([0-9]+(?:\.[0-9]+)?)\s*K\b",
        details,
        re.IGNORECASE,
    )

    if match:
        try:
            return float(match.group(1))
        except ValueError:
            return None

    return None


#.....get_pdbx_ph_range_from_mmcif_or_details...
def get_pdbx_ph_range_from_mmcif_or_details(block):
    """
    Return pH range from mmCIF if present, otherwise extract from pdbx_details.
    Output format: 'low-high' (string), or None if not found.
    """

    # 1️⃣ Try mmCIF field first
    ph_range = block.find_value("_exptl_crystal_grow.pdbx_pH_range")
    if ph_range not in (None, "", ".", "?", "NA"):
        return ph_range.strip()

    # 2️⃣ Fallback: extract from pdbx_details
    details = block.find_value("_exptl_crystal_grow.pdbx_details")
    if not details:
        return None

    details = str(details)

    # Match patterns like:
    # pH 6.5-7.5
    # pH range 6.0 – 7.0
    # pH: 6.8 to 7.2
    match = re.search(
        r"\bpH\b.*?([0-9]+(?:\.[0-9]+)?)\s*(?:-|–|to)\s*([0-9]+(?:\.[0-9]+)?)",
        details,
        re.IGNORECASE,
    )

    if match:
        low = match.group(1)
        high = match.group(2)
        return f"{low}-{high}"

    return None


def extract_mmcif_info(pdb_id, score):
    pdb_id = pdb_id.upper()
    pubmed = get_cached_pubmed_id(pdb_id)
    
    mmcif_content = get_cached_mmcif(pdb_id)
    if not mmcif_content:
        return None
    
    try:
        doc = gemmi.cif.read_string(mmcif_content)
        block = doc.sole_block()
    except Exception as e:
        print(f"⚠ Failed to parse mmCIF for {pdb_id}: {e}")
        return None
    
    info = {
        "PDB_ID": block.find_value("_entry.id"),
        "score": score,
        "Resolution": block.find_value("_refine.ls_d_res_high"),
        "pubmed_id": pubmed,
        "apparatus": block.find_value("_exptl_crystal_grow.apparatus"),
        "atmosphere": block.find_value("_exptl_crystal_grow.atmosphere"),
        "crystal_id ": block.find_value("_exptl_crystal_grow.crystal_id"),
        "details": block.find_value("_exptl_crystal_grow.details"),
        "method": get_method_from_mmcif_or_details(block),
        "method_ref": block.find_value("_exptl_crystal_grow.method_ref"),
        "pH": get_ph_from_mmcif_or_details(block),
        "pressure": block.find_value("_exptl_crystal_grow.pressure"),
        "pressure_esd": block.find_value("_exptl_crystal_grow.pressure_esd"),
        "seeding": block.find_value("_exptl_crystal_grow.seeding"),
        "seeding_ref": block.find_value("_exptl_crystal_grow.seeding_ref"),
        "temp": get_temperature_from_mmcif_or_details(block),
        "temp_esd": block.find_value("_exptl_crystal_grow.temp_esd"),
        "temp_details": block.find_value("_exptl_crystal_grow.temp_details"),
        "time": block.find_value("_exptl_crystal_grow.time"),
        "pdbx_details": block.find_value("_exptl_crystal_grow.pdbx_details"),
        "pdbx_pH_range": get_pdbx_ph_range_from_mmcif_or_details(block),
    }
    return info

def filter_experimental_conditions(input_csv, output_csv=None):
    if output_csv is None:
        base, ext = os.path.splitext(input_csv)
        output_csv = f"{base}_filtered{ext}"
    filtered_rows = []
    with open(input_csv, newline="") as f:
        reader = csv.DictReader(f)
        for row in reader:
            if any(row.get(field) not in (None, "", "NA") for field in ["pH", "temp", "method", "pdbx_details"]):
                filtered_rows.append(row)
    with open(output_csv, "w", newline="") as f:
        writer = csv.DictWriter(f, fieldnames=reader.fieldnames)
        writer.writeheader()
        writer.writerows(filtered_rows)
    print(f"✔ Filtered CSV saved to: {os.path.abspath(output_csv)}")
    return output_csv

def run_pdb_search(sequence, output_csv="pdb_mmcif_extracted.csv", keep_all=False, max_workers=6):
    """
    Run a PDB search for a given protein sequence and extract mmCIF experimental data.
    sequence: str, protein sequence
    keep_all: if True, return full CSV without filtering
    max_workers: number of parallel workers for fetching PDB data (default: 6)
    """

    # Clean sequence
    seq = sequence.strip().upper()
    seq = re.sub(r"\s+", "", seq)
    if not seq:
        raise ValueError("Sequence cannot be empty")
    if not re.fullmatch(r"[ACDEFGHIKLMNPQRSTVWY]+", seq):
        raise ValueError("Sequence contains invalid characters. Use only standard amino acids.")

    # ---------------------------
    # PDB query
    # ---------------------------
    query = {
        "query": {
            "type": "group",
            "logical_operator": "and",
            "nodes": [
                {
                    "type": "terminal",
                    "service": "sequence",
                    "parameters": {
                        "evalue_cutoff": 1e-5,
                        "identity_cutoff": 0.3,
                        "sequence_type": "protein",
                        "value": seq
                    }
                },
                {
                    "type": "terminal",
                    "service": "text",
                    "parameters": {
                        "attribute": "exptl.method",
                        "operator": "exact_match",
                        "value": "X-RAY DIFFRACTION"
                    }
                }
            ]
        },
        "request_options": {"paginate": {"start": 0, "rows": 10000}},
        "return_type": "entry"
    }

    result = requests.post("https://search.rcsb.org/rcsbsearch/v2/query", json=query)
    result.raise_for_status()
    search_data = result.json()
    if search_data["total_count"] == 0:
        raise ValueError("No matching PDB entries found")
    pdb_hits = {hit["identifier"]: hit.get("score") for hit in search_data["result_set"]}

    abs_path = os.path.abspath(output_csv)
    print(f"▶ Writing CSV to: {abs_path}")
    print(f"▶ Found {len(pdb_hits)} PDB entries. Fetching data in parallel ({max_workers} workers)...")

    with open(abs_path, "w", newline="") as f:
        writer = csv.DictWriter(
            f,
            fieldnames=[
                "PDB_ID", "score", "Resolution", "pubmed_id", "apparatus", "atmosphere",
                "crystal_id ", "details", "method", "method_ref", "pH", "pressure",
                "pressure_esd", "seeding", "seeding_ref", "temp", "temp_esd",
                "temp_details", "time", "pdbx_details", "pdbx_pH_range"
            ]
        )
        writer.writeheader()
        
        # Use ThreadPoolExecutor for parallel requests
        with ThreadPoolExecutor(max_workers=max_workers) as executor:
            # Submit all tasks
            future_to_pdb = {
                executor.submit(extract_mmcif_info, pdb_id, score): pdb_id
                for pdb_id, score in pdb_hits.items()
            }
            
            # Process results as they complete
            completed = 0
            for future in as_completed(future_to_pdb):
                completed += 1
                pdb_id = future_to_pdb[future]
                try:
                    row = future.result()
                    if row:  # Only write if extraction was successful
                        writer.writerow(row)
                    if completed % 10 == 0:
                        print(f"  ✓ Processed {completed}/{len(pdb_hits)}")
                except Exception as e:
                    print(f"  ✗ Failed for {pdb_id}: {e}")
            
            print(f"✔ Completed fetching all {len(pdb_hits)} entries")

    if keep_all:
        return abs_path  # return full CSV
    else:
        filtered_csv = filter_experimental_conditions(abs_path)
        return filtered_csv
