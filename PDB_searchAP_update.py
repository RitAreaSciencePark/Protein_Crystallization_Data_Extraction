import requests
import json
import gemmi
import csv
import os


ROW_COUNT = 20

#this query is to navigate the PDB database and find the proteins base on the sequences determined by X-ray diffraction ,methods
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
            "identity_cutoff": 0.9,
            "sequence_type": "protein",
            "value": "MTEYKLVVVGAGGVGKSALTIQLIQNHFVDEYDPTIEDSYRKQVVIDGETCLLDILDTAGQEEYSAMRDQYMRTGEGFLCVFAINNTKSFEDIHQYREQIKRVKDSDDVPMVLVGNKCDLPARTVETRQAQDLARSYGIPYIETSAKTRQGVEDAFYTLVREIRQHKLRKLNPPDESGPGCMNCKCVIS"
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
print("We get a lot of data (and still, is not everything): ")
for labels in entry_data:
  print("\t - "+ labels + ": " + str(entry_data[labels])[:50] + ".......")
print('\n Let\'s focus on the last bit of information:')

#print all external references (e.g., doi, pubmed, etc.)
print(entry_data.get('rcsb_external_references'))

# Loop over all results and print some information about the experimental conditions
search_data = result.json()

if search_data["total_count"] == 0:
    raise ValueError("No matching PDB entries found")

print("\nExtracting experimental conditions:\n")

pdb_ids = [hit["identifier"] for hit in search_data["result_set"]]

def extract_mmcif_info(pdb_id):
    pdb_id = pdb_id.upper()
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
        "PubMed_ID": block.find_value("_citation.pdbx_database_id_PubMed")}

    # ---- Crystallization method (may be multiple) ----
    cryst_methods = []
    for row in block.find_loop("_exptl_crystal_grow"):
        cryst_methods.append(row[0])
    info["apparatus"] = "; ".join(cryst_methods) if cryst_methods else "NA"
    info["atmosphere"] = "; ".join(cryst_methods) if cryst_methods else "NA"
    info["details"] = "; ".join(cryst_methods) if cryst_methods else "NA"
    info["method.refined"] = "; ".join(cryst_methods) if cryst_methods else "NA"
    info["pressure"] = "; ".join(cryst_methods) if cryst_methods else "NA"
    info["pressure_esd"] = "; ".join(cryst_methods) if cryst_methods else "NA"
    info["seeding"] = "; ".join(cryst_methods) if cryst_methods else "NA"
    info["temp_details"] = "; ".join(cryst_methods) if cryst_methods else "NA"
    info["temp_esd"] = "; ".join(cryst_methods) if cryst_methods else "NA"
    info["pdbx_pH_range"] = "; ".join(cryst_methods) if cryst_methods else "NA"
    info["time"] = "; ".join(cryst_methods) if cryst_methods else "NA"
    info["pdbx_details"] = "; ".join(cryst_methods) if cryst_methods else "NA"

    return info

output_csv = "pdb_mmcif_extracted_data.csv"

with open(output_csv, "w", newline="") as f:
    writer = csv.DictWriter(
        f,
        fieldnames=["PDB_ID", "Resolution", "PubMed_ID", "apparatus", "atmosphere", "details", "method.refined", "pressure", "pressure_esd", "seeding", 
    "seeding.ref", "temp", "temp_details", "temp_esd", "time", "pdbx_details", "pdbx_pH_range", "chains"])
    writer.writeheader()

    for pdb_id in pdb_ids:
        try:
            data = extract_mmcif_info(pdb_id)
            writer.writerow(data)
        except Exception as e:
            print(f"Failed for {pdb_id}: {e}")

print(f"\nCSV file saved as: {os.path.abspath(output_csv)}")
