import requests
import csv
import gemmi

RCSB_SEARCH_URL = "https://search.rcsb.org/rcsbsearch/v2/query"
#RCSB_CIF_URL = "https://files.rcsb.org/view/{}.cif"


# ---------- FASTA READER ----------
def read_fasta(fasta_file):
    sequences = {}
    header = None
    seq_parts = []

    with open(fasta_file) as f:
        for line in f:
            line = line.strip()
            if line.startswith(">"):
                if header:
                    sequences[header] = "".join(seq_parts)
                header = line[1:].split()[0]
                seq_parts = []
            else:
                seq_parts.append(line)
        if header:
            sequences[header] = "".join(seq_parts)

    return sequences

# ---------- PDB SEARCH ----------
def search_pdb(sequence, identity_cutoff=0.5, rows=10):

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
            "paginate": {
                "start": 0,
                "rows": rows
            }
        },
        "return_type": "entry"
    }
    query["request_options"]["paginate"]["rows"] = rows

    r = requests.post(RCSB_SEARCH_URL, json=query, timeout=10)
    r.raise_for_status()

    return [hit["identifier"] for hit in r.json()["result_set"]]

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

    cif_file = f"{pdb_id}.cif"
    with open(cif_file, "wb") as f:
        f.write(response.content)

    doc = gemmi.cif.read(cif_file)
    block = doc.sole_block()

    return {
        "PDB_ID": block.find_value("_entry.id"),
            "Resolution": block.find_value("_refine.ls_d_res_high"),
            "pubmed_id": pubmed,
            "apparatus": block.find_value("_exptl_crystal_grow.apparatus"),
            "atmosphere": block.find_value("_exptl_crystal_grow.atmosphere"),
            "crystal_id ": block.find_value("_exptl_crystal_grow.crystal_id"),
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


# ---------- MAIN PIPELINE ----------
def fasta_to_pdb_csv(fasta_file, output_csv):
    sequences = read_fasta(fasta_file)

    with open(output_csv, "w", newline="") as f:

        raw_writer = csv.writer(f)

        fieldnames=["PDB_ID", "Resolution", "pubmed_id", "apparatus", "atmosphere", "crystal_id ", "details", 
            "method", "method_ref", "pH", "pressure", "pressure_esd", "seeding", "seeding_ref", "temp", "temp_esd", 
            "temp_details", "time", "pdbx_details", "pdbx_pH_range"]

        writer = csv.DictWriter(f, fieldnames=fieldnames)
        writer.writeheader()


        for query_id, sequence in sequences.items():
            print(f"Searching PDB for {query_id}...")

            #...Write the sequence line no column
            raw_writer.writerow([query_id, sequence])

            pdb_ids = search_pdb(sequence)

            for pdb_id in pdb_ids:
                try:
                    cryst = extract_crystallization(pdb_id)
                    writer.writerow({
                        "PDB_ID": pdb_id,
                        **cryst
                    })
                except Exception as e:
                    print(f"Failed for {pdb_id}: {e}")


# ---------- RUN ----------
if __name__ == "__main__":
    fasta_to_pdb_csv("cleaned_KaiB_MSA.fasta", "cleaned_KaiB_MSA_crystallization_results.csv")
