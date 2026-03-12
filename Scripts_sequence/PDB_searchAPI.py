import requests
import json
import gemmi
import csv
import os
import re
import sys
from concurrent.futures import ThreadPoolExecutor, as_completed
from Bio.PDB.MMCIF2Dict import MMCIF2Dict

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

def fetch_ligands(pdb_id):
    """
    Fetch ligand IDs for a PDB entry.

    1️⃣ Try RCSB GraphQL
    2️⃣ If none found → fallback to mmCIF parsing
    """

    pdb_id = pdb_id.upper()

    # ---------- GraphQL ----------
    url = "https://data.rcsb.org/graphql"

    query = """
    query getLigands($id: String!) {
      entry(entry_id: $id) {
        nonpolymer_entities {
          rcsb_nonpolymer_entity_container_identifiers {
            nonpolymer_comp_id
          }
        }
      }
    }
    """

    try:
        response = requests.post(
            url,
            json={"query": query, "variables": {"id": pdb_id}},
            timeout=15
        )

        response.raise_for_status()
        data = response.json()

        entry = data.get("data", {}).get("entry", {})
        entities = entry.get("nonpolymer_entities") or []

        ligands = [
            e["rcsb_nonpolymer_entity_container_identifiers"]["nonpolymer_comp_id"]
            for e in entities
            if e.get("rcsb_nonpolymer_entity_container_identifiers", {}).get("nonpolymer_comp_id")
        ]

        if ligands:
            return ", ".join(list(dict.fromkeys(ligands)))

    except Exception as e:
        print(f"⚠ GraphQL ligand error for {pdb_id}: {e}")

    # ---------- mmCIF fallback ----------
    try:

        cif_path = os.path.join(CACHE_DIR, f"{pdb_id.lower()}.cif")

        if not os.path.exists(cif_path):
            mmcif_content = get_cached_mmcif(pdb_id)
            if not mmcif_content:
                return None

        doc = gemmi.cif.read_file(cif_path)
        block = doc.sole_block()

        comp_ids = block.find_values("_pdbx_entity_instance_feature.comp_id")

        comp_ids = [c for c in comp_ids if c not in ("?", ".", "", None)]

        if comp_ids:
            return ", ".join(list(dict.fromkeys(comp_ids)))

    except Exception as e:
        print(f"⚠ mmCIF ligand fallback error for {pdb_id}: {e}")

    return None

def extract_mmcif_info(pdb_id, score):
    pdb_id = pdb_id.upper()

    # 1️⃣ PubMed
    pubmed = get_cached_pubmed_id(pdb_id)

    # 2️⃣ Download mmCIF
    mmcif_content = get_cached_mmcif(pdb_id)
    if not mmcif_content:
        return None

    # 3️⃣ Parse mmCIF with gemmi
    try:
        doc = gemmi.cif.read_string(mmcif_content)
        block = doc.sole_block()
    except Exception as e:
        print(f"⚠ Failed to parse mmCIF for {pdb_id}: {e}")
        return None

    # ---- Extract assembly info robustly ----
    assembly_ids = block.find_values("_pdbx_struct_assembly.id")
    oligomeric_values = block.find_values("_pdbx_struct_assembly.oligomeric_details")

    if not oligomeric_values:
        assembly_detail = None
    else:
        first_detail = oligomeric_values[0]
        if len(assembly_ids) > 1:
            first_detail += "(*)"
        assembly_detail = first_detail

    # 5️⃣ Collect info
    info = {
        "PDB_ID": block.find_value("_entry.id"),
        "score": score,
        "pubmed_id": pubmed,
        "Assembly": assembly_detail,
        "Resolution": block.find_value("_refine.ls_d_res_high"),
        "pdbx_details": block.find_value("_exptl_crystal_grow.pdbx_details"),
        "method": get_method_from_mmcif_or_details(block),
        "pH": get_ph_from_mmcif_or_details(block),
        "temp": get_temperature_from_mmcif_or_details(block),
        "pdbx_pH_range": get_pdbx_ph_range_from_mmcif_or_details(block),
        "ligands": fetch_ligands(pdb_id)
    }

    return info

def filter_experimental_conditions(input_csv, output_csv=None):

   
    if output_csv is None:
        base, ext = os.path.splitext(input_csv)
        output_csv = f"{base}_filtered{ext}"
    filtered_rows = []

    if os.path.exists(input_csv):
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
    else:
       print(f" CSV not found: {input_csv}, skipping processing.")
       sys.exit(1)

def detect_seq_type(sequence):

    """Automatically detect whether a sequence is protein, DNA, or RNA."""
    sequence = sequence.upper()
    dna_letters = set("ACGT")
    rna_letters = set("ACGU")
    protein_letters = set("ACDEFGHIKLMNPQRSTVWY")  # standard amino acids

    seq_set = set(sequence)

    if seq_set <= protein_letters:
        return "protein"
    elif seq_set <= dna_letters:
        return "dna"
    elif seq_set <= rna_letters:
        return "rna"
    else:
        # Mixed or unknown letters → default to protein
        return "protein"


def search_pdb_by_sequence(sequence, output_csv="pdb_mmcif_extracted.csv", keep_all=False, max_workers=6):
    """Search PDB by sequence and automatically detect sequence type."""
    
    # Automatically detect the sequence type
    seq_type = detect_seq_type(sequence)

    target_map = {
        "protein": "pdb_protein_sequence",
        "dna": "pdb_dna_sequence",
        "rna": "pdb_rna_sequence"
    }

    target = target_map[seq_type]

    query = {
        "query": {
            "type": "group",
            "logical_operator": "and",
            "nodes": [

                # Sequence similarity search
                {
                    "type": "terminal",
                    "service": "sequence",
                    "parameters": {
                        "target": target,
                        "value": sequence,
                        "identity_cutoff": 0.3,
                        "evalue_cutoff": 1e-3
                    }
                },

                # Filter: only X-ray diffraction structures
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

    # The rest of your function can remain unchanged

    headers = {"User-Agent": "PDB-sequence-search-script/1.0"}

    try:
        result = requests.post(
            "https://search.rcsb.org/rcsbsearch/v2/query",
            json=query,
            headers=headers,
            timeout=30
        )
        result.raise_for_status()
    except requests.exceptions.RequestException as e:
        print(f"⚠ PDB search request failed: {e}")
        return []

    # Check for empty response
    if not result.text.strip():
        print("⚠ PDB API returned empty response")
        return []

    # Parse JSON safely
    try:
        search_data = result.json()
    except ValueError:
        print("⚠ PDB API returned invalid JSON. Raw response snippet:")
        print(result.text[:500])
        return []

    # Safe check for results
    if search_data.get("total_count", 0) == 0:
        print("⚠ No matching PDB entries found")
        return []

    pdb_hits = {hit["identifier"]: hit.get("score") for hit in search_data.get("result_set", [])}

    if not pdb_hits:
        print("❌ No hits for this sequence. Stopping pipeline.")
        return None

    # Always write CSV with headers, even if empty
    fieldnames = [
        "PDB_ID","score","Resolution","pubmed_id","Assembly","method","pH","temp","pdbx_details","pdbx_pH_range","ligands"]
    
    os.makedirs(os.path.dirname(output_csv), exist_ok=True)
    with open(output_csv, "w", newline="") as f:
        writer = csv.DictWriter(f, fieldnames=fieldnames)
        writer.writeheader()

    if not pdb_hits:
        print("⚠ No PDB entries found. Empty CSV created.")
        return output_csv  # CSV exists with headers

    print(f"▶ Found {len(pdb_hits)} PDB entries. Fetching mmCIF data in parallel ({max_workers} workers)...")

    with open(output_csv, "a", newline="") as f:
        writer = csv.DictWriter(f, fieldnames=fieldnames)
        with ThreadPoolExecutor(max_workers=max_workers) as executor:
            future_to_pdb = {executor.submit(extract_mmcif_info, pdb_id, score): pdb_id for pdb_id, score in pdb_hits.items()}
            completed = 0
            for future in as_completed(future_to_pdb):
                completed += 1
                pdb_id = future_to_pdb[future]
                try:
                    row = future.result()
                    if row:
                        writer.writerow(row)
                    if completed % 10 == 0:
                        print(f"  ✓ Processed {completed}/{len(pdb_hits)}")
                except Exception as e:
                    print(f"  ✗ Failed for {pdb_id}: {e}")

    print(f"✔ CSV written to: {os.path.abspath(output_csv)}")
    return output_csv
