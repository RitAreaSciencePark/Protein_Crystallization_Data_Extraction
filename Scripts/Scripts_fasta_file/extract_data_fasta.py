import requests
import csv
import gemmi
import time
import random
import argparse



RCSB_SEARCH_URL = "https://search.rcsb.org/rcsbsearch/v2/query"
#RCSB_CIF_URL = "https://files.rcsb.org/view/{}.cif"

def read_fasta(fasta_file):
    fasta_dict = {}
    current_id = None
    seq_lines = []

    with open(fasta_file) as f:
        for line in f:
            line = line.strip()
            if not line:
                continue
            if line.startswith(">"):
                if current_id:
                    fasta_dict[current_id] = "".join(seq_lines)
                     # Save previous sequence
                current_id = line[1:].split()[0]  # Take entire header or first word as ID
                seq_lines = []
            else:
                seq_lines.append(line)

        if current_id:
            fasta_dict[current_id] = "".join(seq_lines)

    return fasta_dict

def extract_pubmed_id(pdb_id):
    entry_data = requests.get(f"https://data.rcsb.org/rest/v1/core/entry/{pdb_id}").json()

    pubmed = entry_data.get("rcsb_primary_citation", {}).get("pdbx_database_id_pub_med", "NA")
    return pubmed

# ---------- mmCIF EXTRACTION ----------
def extract_crystallization(pdb_id):
    #pdb_id = pdb_id.upper()
    pubmed = extract_pubmed_id(pdb_id)

    response = requests.get(f"https://files.rcsb.org/view/{pdb_id}.cif")
    response.raise_for_status()

    doc = gemmi.cif.read_string(response.text)
    block = doc.sole_block()

    return {
            "Resolution": block.find_value("_refine.ls_d_res_high"),
            "pubmed_id": pubmed,
            "apparatus": block.find_value("_exptl_crystal_grow.apparatus"),
            "atmosphere": block.find_value("_exptl_crystal_grow.atmosphere"),
            "crystal_id": block.find_value("_exptl_crystal_grow.crystal_id"),
            "details": block.find_value("_exptl_crystal_grow.details"),
            "method": block.find_value("_exptl_crystal_grow.method"),
            "method_ref": block.find_value("_exptl_crystal_grow.method_ref"),
            "pH": block.find_value("_exptl_crystal_grow.pH"),
            "pressure": block.find_value("_exptl_crystal_grow.pressure"),
            "pressure_esd": block.find_value("_exptl_crystal_grow.pressure_esd"),
            "seeding": block.find_value("_exptl_crystal_grow.seeding"),
            "seeding_ref": block.find_value("_exptl_crystal_grow.seeding_ref"),
            "temp": block.find_value("_exptl_crystal_grow.temp"),
            "temp_esd": block.find_value("_exptl_crystal_grow.temp_esd"),
            "temp_details": block.find_value("_exptl_crystal_grow.temp_details"),
            "time": block.find_value("_exptl_crystal_grow.time"),
            "pdbx_details": block.find_value("_exptl_crystal_grow.pdbx_details"),
            "pdbx_pH_range": block.find_value("_exptl_crystal_grow.pdbx_pH_range")}


# --- Search PDB and save all results to one CSV ---
def fasta_to_pdb_csv(fasta_file, output_csv, identity_cutoff=0.2, rows=1000, delay=1):
    sequences = read_fasta(fasta_file)

    fieldnames = [
        "query_id", "pdb_id", "score", "Resolution", "pubmed_id", "apparatus", "atmosphere",
        "crystal_id", "details", "method", "method_ref", "pH", "pressure", "pressure_esd",
        "seeding", "seeding_ref", "temp", "temp_esd", "temp_details", "time", "pdbx_details", "pdbx_pH_range"
    ]

    with open(output_csv, "w", newline="") as f:
        writer = csv.DictWriter(f, fieldnames=fieldnames)
        writer.writeheader()

        for query_id, sequence in sequences.items():
            query = {
                "query": {
                    "type": "group",
                    "logical_operator": "and",
                    "nodes": [
                        {
                            "type": "terminal",
                            "service": "sequence",
                            "parameters": {
                                "sequence_type": "protein",
                                "value": sequence,
                                "identity_cutoff": identity_cutoff,
                                "evalue_cutoff": 1
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
                "request_options": {
                    "paginate": {"start": 0, "rows": rows}
                },
                "return_type": "entry"
            }

            try:
                response = requests.post(RCSB_SEARCH_URL, json=query, timeout=30)
                response.raise_for_status()
                if not response.text.strip():
                    print(f"[WARN] {query_id}: Empty response")
                    continue
                data = response.json()
            except Exception as e:
                print(f"[ERROR] {query_id}: {e}")
                continue

            total_hits = data.get("total_count", 0)
            if total_hits == 0:
                print(f"{query_id}: 0 hits, skipping")
                continue

            print(f"{query_id}: {total_hits} hits")

            for entry in data.get("result_set", []):
                pdb_id = entry.get("identifier")
                score = entry.get("score")
                crystallization = extract_crystallization(pdb_id)

                row = {"query_id": query_id, "pdb_id": pdb_id, "score": score}
                row.update(crystallization)
                writer.writerow(row)
                f.flush()

           # time.sleep(delay)

    print(f"PDB search results saved to {output_csv}")


