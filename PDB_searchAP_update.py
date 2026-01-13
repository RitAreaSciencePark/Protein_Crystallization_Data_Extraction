import requests
import json
import gemmi
import csv
import os
import re

def get_sequence_from_user():
    seq = input(
        "Enter protein sequence (single-letter amino acid code, no FASTA header):\n"
    ).strip().upper()

    # Remove whitespace and line breaks
    seq = re.sub(r"\s+", "", seq)

    if not seq:
        raise ValueError("Sequence cannot be empty")

    # Allowed amino acids for RCSB sequence search
    if not re.fullmatch(r"[ACDEFGHIKLMNPQRSTVWY]+", seq):
        raise ValueError(
            "Sequence contains invalid characters. "
            "Use only standard amino acids (ACDEFGHIKLMNPQRSTVWY)."
        )

    return seq


def run_pdb_search(identity_cutoff=0.9, row_count=20):
    """
    Run a PDB search for a given sequence and extract mmCIF experimental data.
    """
    # ✅ Ask for sequence HERE
    seq = get_sequence_from_user()
    # ---------- OUTPUT FILE ----------
    output_csv = "pdb_mmcif_extracted_data.csv"
    ROW_COUNT = 20


#this query is to navigate the PDB database and find the proteins base on the sequences determined by X-ray diffraction ,methods
    query = {
    "query": {
        "type": "group",
        "logical_operator": "and",
        "nodes": [
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
            "rows": 10
        }
        },
        "return_type": "entry"
    }

    # Load the query string/JSON into a PYTHON dictionary
    #query = json.loads(query)
    query["request_options"]["paginate"]["rows"] = ROW_COUNT

    # Post our query to the system
    result = requests.post("https://search.rcsb.org/rcsbsearch/v2/query",json=query)

    # If you want to know the total number of pages (e.g., to do an automated extraction)
    result.raise_for_status()  # stops if request failed
    PAGE_MAX = result.json()["total_count"]//ROW_COUNT + 1
    print("rows: ", result.json()["total_count"], "pages:", PAGE_MAX)

    # Let's select the first result:

    if result.json()["total_count"] > 0:
        print("First result:", result.json()["result_set"][0])
    else:
        print("No matching structures found.")

    print(result)
    print(result.json())

    #Get a PDB ID from the search results
    identifier = result.json()["result_set"][0]["identifier"]

    #Fetch detailed metadata for that PDB entry
    entry_data = requests.get(f'https://data.rcsb.org/rest/v1/core/entry/{identifier}').json()
   

    #Print all top-level data fields (preview)
    #print("We get a lot of data (and still, is not everything): ")
    #for labels in entry_data:
    #print("\t - "+ labels + ": " + str(entry_data[labels])[:50] + ".......")
    #print('\n Let\'s focus on the last bit of information:')

    #print all external references (e.g., doi, pubmed, etc.)
    #print(entry_data.get('rcsb_external_references'))

    # Loop over all results and print some information about the experimental conditions
    search_data = result.json()

    if search_data["total_count"] == 0:
        raise ValueError("No matching PDB entries found")

    print("\nExtracting experimental conditions:\n")

    pdb_ids = [hit["identifier"] for hit in search_data["result_set"]]

    def extract_pubmed_id(pdb_id):
        """
        Fetch PubMed ID for a given PDB entry.
        """
        entry_data = requests.get(
            f"https://data.rcsb.org/rest/v1/core/entry/{pdb_id}"
        ).json()

        pubmed = entry_data.get(
            "rcsb_primary_citation", {}
        ).get(
            "pdbx_database_id_pub_med", "NA"
        )

        return pubmed


    def extract_mmcif_info(pdb_id):
        pdb_id = pdb_id.upper()

        pubmed = extract_pubmed_id(pdb_id)
        cif_file = f"{pdb_id}.cif"

        # ---- Download mmCIF file ----
        response = requests.get(f"https://files.rcsb.org/view/{pdb_id}.cif")
        response.raise_for_status()

        with open(cif_file, "wb") as f:
            f.write(response.content)

        # ---- Read mmCIF ----
        doc = gemmi.cif.read(cif_file)
        block = doc.sole_block()     

        info = {
            "PDB_ID": block.find_value("_entry.id"),
            "Resolution": block.find_value("_refine.ls_d_res_high"),
            "pubmed_id": pubmed,
            "apparatus": block.find_value("_exptl_crystal_grow.apparatus"),
            "atmosphere": block.find_value("_exptl_crystal_grow.atmosphere"),
            "crystal_id ": block.find_value("_exptl_crystal_grow.crystal_id"),
            "details": block.find_value("_exptl_crystal_grow.details"),
            "method": block.find_value("_exptl_crystal_grow.methods"),
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
            "pdbx_pH_range": block.find_value("_exptl_crystal_grow.pdbx_pH_range"),}

        return info


    output_csv = "pdb_mmcif_extracted.csv"

    with open(output_csv, "w", newline="") as f:
        writer = csv.DictWriter(
            f,
            fieldnames=["PDB_ID", "Resolution", "pubmed_id", "apparatus", "atmosphere", "crystal_id ", "details", 
            "method", "method_ref", "pH", "pressure", "pressure_esd", "seeding", "seeding_ref", "temp", "temp_esd", 
            "temp_details", "time", "pdbx_details", "pdbx_pH_range"])
        writer.writeheader()

        for pdb_id in pdb_ids:
            try:
                data = extract_mmcif_info(pdb_id)
                writer.writerow(data)
            except Exception as e:
                print(f"Failed for {pdb_id}: {e}")

    print(f"\nCSV file saved as: {os.path.abspath(output_csv)}")

if __name__ == "__main__":
    run_pdb_search()
