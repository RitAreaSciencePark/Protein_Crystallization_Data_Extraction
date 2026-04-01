import os
import csv
import re
import requests
import gemmi
import pandas as pd
from concurrent.futures import ThreadPoolExecutor, as_completed
from Bio.Align import PairwiseAligner

# ------------------ CACHE SETUP ------------------
CACHE_DIR = os.path.join(os.path.dirname(os.path.abspath(__file__)), ".pdb_cache")
os.makedirs(CACHE_DIR, exist_ok=True)

# ------------------ HELPER FUNCTIONS ------------------
def get_cached_pubmed_id(pdb_id):
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
        with open(cache_file, "w") as f:
            json.dump({"pubmed_id": pubmed}, f)
        return pubmed
    except Exception as e:
        print(f"⚠ Warning fetching PubMed ID for {pdb_id}: {e}")
        return "NA"

def get_cached_mmcif(pdb_id):
    cache_file = os.path.join(CACHE_DIR, f"{pdb_id.lower()}.cif")
    if os.path.exists(cache_file):
        with open(cache_file, "r") as f:
            return f.read()
    try:
        url = f"https://files.rcsb.org/view/{pdb_id}.cif"
        response = requests.get(url, timeout=10)
        response.raise_for_status()
        content = response.text
        with open(cache_file, "w") as f:
            f.write(content)
        return content
    except Exception as e:
        print(f"⚠ Warning fetching mmCIF for {pdb_id}: {e}")
        return None

def detect_seq_type(sequence):
    sequence = sequence.upper()
    dna_letters = set("ACGT")
    rna_letters = set("ACGU")
    protein_letters = set("ACDEFGHIKLMNPQRSTVWY")
    seq_set = set(sequence)
    if seq_set <= protein_letters:
        return "protein"
    elif seq_set <= dna_letters:
        return "dna"
    elif seq_set <= rna_letters:
        return "rna"
    else:
        return "protein"

def compute_identity(seq1, seq2):
    """Compute sequence identity using global alignment with PairwiseAligner"""
    aligner = PairwiseAligner()
    aligner.mode = "global"
    aligner.match_score = 1
    aligner.mismatch_score = 0
    aligner.open_gap_score = 0
    aligner.extend_gap_score = 0
    score = aligner.score(seq1, seq2)
    alignment_length = max(len(seq1), len(seq2))
    identity = (score / alignment_length) * 100
    return identity

def fetch_pdb_sequences(pdb_id, seq_type="protein"):
    pdb_id = pdb_id.upper()
    mmcif_content = get_cached_mmcif(pdb_id)
    if not mmcif_content:
        return []
    try:
        doc = gemmi.cif.read_string(mmcif_content)
        block = doc.sole_block()
    except Exception as e:
        print(f"⚠ Failed to parse mmCIF for sequences {pdb_id}: {e}")
        return []

    sequences = []
    if seq_type == "protein":
        seq_category = "_entity_poly.pdbx_seq_one_letter_code_can"
    elif seq_type in ("dna", "rna"):
        seq_category = "_entity_poly.pdbx_seq_one_letter_code"
    else:
        return []

    seqs = block.find_values(seq_category)
    for s in seqs:
        if s and s not in ("?", ".", ""):
            sequences.append(s.replace("\n", "").replace(" ", ""))
    return sequences

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

def get_temperature_from_mmcif_or_details(block):
    temp = block.find_value("_exptl_crystal_grow.temp")
    if temp not in (None, "", ".", "?", "NA"):
        try:
            return float(temp)
        except ValueError:
            pass
    details = block.find_value("_exptl_crystal_grow.pdbx_details")
    if not details:
        return None
    match = re.search(r"\b([0-9]+(?:\.[0-9]+)?)\s*K\b", details, re.IGNORECASE)
    if match:
        try:
            return float(match.group(1))
        except ValueError:
            return None
    return None

def get_pdbx_ph_range_from_mmcif_or_details(block):
    ph_range = block.find_value("_exptl_crystal_grow.pdbx_pH_range")
    if ph_range not in (None, "", ".", "?", "NA"):
        return ph_range.strip()
    details = block.find_value("_exptl_crystal_grow.pdbx_details")
    if not details:
        return None
    match = re.search(r"\bpH\b.*?([0-9]+(?:\.[0-9]+)?)\s*(?:-|–|to)\s*([0-9]+(?:\.[0-9]+)?)", str(details), re.IGNORECASE)
    if match:
        return f"{match.group(1)}-{match.group(2)}"
    return None

def fetch_ligands(pdb_id):
    pdb_id = pdb_id.upper()
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
        print(f"⚠ Ligand fetch error for {pdb_id}: {e}")
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

    assembly_ids = block.find_values("_pdbx_struct_assembly.id")
    oligomeric_values = block.find_values("_pdbx_struct_assembly.oligomeric_details")
    if not oligomeric_values:
        assembly_detail = None
    else:
        first_detail = oligomeric_values[0]
        if len(assembly_ids) > 1:
            first_detail += "(*)"
        assembly_detail = first_detail

    info = {
        "PDB_ID": block.find_value("_entry.id"),
        "Score": score,  # ← use computed max identity
        "Pubmed_id": pubmed,
        "Assembly": assembly_detail,
        "Resolution": block.find_value("_refine.ls_d_res_high"),
        "pdbx_details": block.find_value("_exptl_crystal_grow.pdbx_details"),
        "Method": get_method_from_mmcif_or_details(block),
        "pH": get_ph_from_mmcif_or_details(block),
        "Temp": get_temperature_from_mmcif_or_details(block),
        "pdbx_pH_range": get_pdbx_ph_range_from_mmcif_or_details(block),
        "Ligands": fetch_ligands(pdb_id)
    }
    return info

# ------------------ MAIN SEARCH FUNCTION ------------------
def search_pdb_by_sequence(sequence, output_csv="pdb_mmcif_extracted.csv", max_workers=6):
    seq_type = detect_seq_type(sequence)
    target_map = {"protein":"pdb_protein_sequence","dna":"pdb_dna_sequence","rna":"pdb_rna_sequence"}
    target = target_map[seq_type]

    query = {
        "query": {
            "type": "group",
            "logical_operator": "and",
            "nodes": [
                {"type":"terminal","service":"sequence","parameters":{"target":target,"value":sequence,"identity_cutoff":0.5,"evalue_cutoff":1e-5}},
                {"type":"terminal","service":"text","parameters":{"attribute":"exptl.method","operator":"exact_match","value":"X-RAY DIFFRACTION"}}
            ]
        },
        "request_options":{"paginate":{"start":0,"rows":10000}},
        "return_type":"entry"
    }

    headers = {"User-Agent": "PDB-sequence-search-script/1.0"}
    try:
        result = requests.post("https://search.rcsb.org/rcsbsearch/v2/query", json=query, headers=headers, timeout=30)
        result.raise_for_status()
        search_data = result.json()
    except Exception as e:
        print(f"⚠ PDB search request failed or invalid response: {e}")
        return []

    if search_data.get("total_count",0) == 0:
        print("⚠ No matching PDB entries found")
        return []

    pdb_hits = {hit["identifier"]: hit.get("score", 0) for hit in search_data.get("result_set", [])}
    if not pdb_hits:
        print("❌ No hits for this sequence. Stopping pipeline.")
        return []

    # ---------- Compute max sequence identity ----------
    for pdb_id in list(pdb_hits.keys()):
        pdb_seqs = fetch_pdb_sequences(pdb_id, seq_type=seq_type)
        if pdb_seqs:
            max_identity = max(compute_identity(sequence.upper(), pdb_seq.upper()) for pdb_seq in pdb_seqs)
            pdb_hits[pdb_id] = max_identity
        else:
            pdb_hits[pdb_id] = 0.0

    # ---------- Sort by Score descending ----------
    sorted_pdb_hits = dict(sorted(pdb_hits.items(), key=lambda item: item[1], reverse=True))

    # ---------- Write CSV ----------
    fieldnames = ["PDB_ID","Score","Resolution","Pubmed_id","Assembly","Method","pH","Temp","pdbx_details","pdbx_pH_range","Ligands"]
    os.makedirs(os.path.dirname(output_csv), exist_ok=True)
    with open(output_csv, "w", newline="") as f:
        writer = csv.DictWriter(f, fieldnames=fieldnames)
        writer.writeheader()

    print(f"▶ Found {len(sorted_pdb_hits)} PDB entries. Fetching mmCIF data in parallel ({max_workers} workers)...")
    with open(output_csv, "a", newline="") as f:
        writer = csv.DictWriter(f, fieldnames=fieldnames)
        with ThreadPoolExecutor(max_workers=max_workers) as executor:
            future_to_pdb = {executor.submit(extract_mmcif_info, pdb_id, score): pdb_id 
                             for pdb_id, score in sorted_pdb_hits.items()}
            completed = 0
            for future in as_completed(future_to_pdb):
                completed += 1
                pdb_id = future_to_pdb[future]
                try:
                    row = future.result()
                    if row:
                        writer.writerow(row)
                    if completed % 10 == 0:
                        print(f"  ✓ Processed {completed}/{len(sorted_pdb_hits)}")
                except Exception as e:
                    print(f"  ✗ Failed for {pdb_id}: {e}")

    print(f"✔ CSV written to: {os.path.abspath(output_csv)}")
    return output_csv