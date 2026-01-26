import requests
import gemmi
import csv
import os
import re

def get_sequence_from_user():
    seq = input(
        "Enter protein sequence (single-letter amino acid code, no FASTA header):\n"
    ).strip().upper()

    seq = re.sub(r"\s+", "", seq)
    if not seq:
        raise ValueError("Sequence cannot be empty")

    if not re.fullmatch(r"[ACDEFGHIKLMNPQRSTVWY]+", seq):
        raise ValueError(
            "Sequence contains invalid characters. "
            "Use only standard amino acids (ACDEFGHIKLMNPQRSTVWY)."
        )

    return seq

def parse_from_pdbx_details(details_text):
    """
    Extract crystallization conditions from pdbx_details free text.
    Returns a dictionary with unique buffers, salts, PEG, protein concentration, etc.
    """
    result = {
        "protein_concentration": None,
        "peg": None,
        "pH": None,
        "pH_range": None,
        "time": None,
        "seeding": None,
        "ligand": None,
        "solvent": None
    }

    if not details_text or details_text in {".", "?"}:
        return result

    text = details_text.replace("\n", " ").lower()

    # --- Protein concentration + protein name ---
    protein_matches = re.findall(
        r"(?:protein(?: mixture)?[:\s]*)?(\d+(?:\.\d+)?\s*(?:mg/ml|g/l|µm|um|mm))\s+([a-z0-9\-\s]+?)(?=,|;|$|\bin\b)",
        text
    )
    proteins = []
    for conc, name in protein_matches:
        # remove buffer/salt keywords
        name_clean = re.sub(
            r"\b(tris|hepes|mes|mops|cacodylate|acetate|phosphate|citrate|nacl|kcl|mgcl2|cacl2|zncl2|ammonium sulfate|lithium sulfate|sodium sulfate|potassium chloride)\b",
            "", name
        ).strip()
        if name_clean:
            proteins.append(f"{conc} {name_clean}")
    if proteins:
        result["protein_concentration"] = ", ".join(proteins)

    # --- PEG with percentage ---
    peg_match = re.search(r"(?:peg|polyethylene glycol)\s*(\d{1,3}(?:\.\d+)?)\s*%?", text)
    if peg_match:
        result["peg"] = f"{peg_match.group(1)}%"

    # --- pH ---
    ph_match = re.search(r"ph\s*[=:]?\s*([0-9]+(?:\.[0-9]+)?)", text)
    if ph_match:
        result["pH"] = ph_match.group(1)

    # --- pH range ---
    ph_range_match = re.search(r"ph\s*[=:]?\s*([0-9]+(?:\.[0-9]+)?)\s*(?:-|–|to)\s*([0-9]+(?:\.[0-9]+)?)", text)
    if ph_range_match:
        result["pH_range"] = f"{ph_range_match.group(1)}-{ph_range_match.group(2)}"

    # --- Time ---
    time_match = re.search(r"(\d+(?:\.\d+)?)\s*(days?|hours?|hrs?|h|weeks?)", text)
    if time_match:
        result["time"] = time_match.group(0)

    # --- Seeding ---
    seeding_match = re.search(r"(?:microseeding|macroseeding|seeding|seeded|unseeded|no seeding)", text)
    if seeding_match:
        result["seeding"] = seeding_match.group(0)

     # --- Ligand ---
    ligand_match = re.search(r"(ligand|inhibitor|substrate|cofactor|with\s+[a-z0-9\-]+)", text)
    if ligand_match:
        result["ligand"] = ligand_match.group(0)

    # --- Solvent ---
    solvent_match = re.search(r"\b(water|h2o|ethanol|isopropanol|mpd|glycerol|dmso)\b", text)
    if solvent_match:
        result["solvent"] = solvent_match.group(1)

    return result


def extract_pubmed_id(pdb_id):
    entry_data = requests.get(f"https://data.rcsb.org/rest/v1/core/entry/{pdb_id}").json()
    return entry_data.get("rcsb_primary_citation", {}).get("pdbx_database_id_pub_med", "NA")

def extract_mmcif_info(pdb_id, score):
    pdb_id = pdb_id.upper()
    pubmed = extract_pubmed_id(pdb_id)

    url = f"https://files.rcsb.org/view/{pdb_id}.cif"
    response = requests.get(url)
    response.raise_for_status()
    doc = gemmi.cif.read_string(response.text)
    block = doc.sole_block()

    pdbx_details = block.find_value("_exptl_crystal_grow.pdbx_details")
    pdbx_ph_range = block.find_value("_exptl_crystal_grow.pdbx_pH_range")

    extracted = parse_from_pdbx_details(pdbx_details)

    info =  info = {
        "PDB_ID": block.find_value("_entry.id"),
        "score": score,
        "Resolution": block.find_value("_refine.ls_d_res_high"),
        "pubmed_id": pubmed,
        "details": pdbx_details,
        "method": block.find_value("_exptl_crystal_grow.method"),
        "pH": extracted["pH"],
        "pH_range": pdbx_ph_range,
        "time": extracted["time"],
        "seeding": extracted["seeding"],
        "protein_concentration": extracted["protein_concentration"],
        "peg": extracted["peg"],
        "ligand": extracted["ligand"],
        "solvent": extracted["solvent"]
    }
    return info

def run_pdb_search(identity_cutoff=0.9, row_count=20):
    seq = get_sequence_from_user()
    output_csv = "pdbx_details_extracted.csv"

    query = {
        "query": {
            "type": "group",
            "logical_operator": "and",
            "nodes": [
                {
                    "type": "terminal",
                    "service": "sequence",
                    "parameters": {
                        "evalue_cutoff": 1,
                        "identity_cutoff": identity_cutoff,
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
        "request_options": {
            "paginate": {"start": 0, "rows": row_count}
        },
        "return_type": "entry"
    }

    result = requests.post("https://search.rcsb.org/rcsbsearch/v2/query", json=query)
    result.raise_for_status()
    search_data = result.json()

    if search_data["total_count"] == 0:
        print("No matching PDB entries found.")
        return

    pdb_hits = {hit["identifier"]: hit.get("score") for hit in search_data["result_set"]}

    with open(output_csv, "w", newline="") as f:
        fieldnames = ["PDB_ID", "score", "Resolution", "pubmed_id", "details",
                "method", "pH", "pH_range", "time", "seeding", "protein_concentration", "peg", "ligand", "solvent"]
        writer = csv.DictWriter(f, fieldnames=fieldnames)
        writer.writeheader()

        for pdb_id, score in pdb_hits.items():
            try:
                row = extract_mmcif_info(pdb_id, score)
                writer.writerow(row)
            except Exception as e:
                print(f"Failed for {pdb_id}: {e}")

    print(f"\nCSV file saved as: {os.path.abspath(output_csv)}")


if __name__ == "__main__":
    run_pdb_search()
