import requests
import json


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

for hit in search_data["result_set"]:
    pdb_id = hit["identifier"]

    # Fetch entry metadata
    entry_data = requests.get(
  f"https://data.rcsb.org/rest/v1/core/entry/{pdb_id}").json()

    # Crystallization conditions (may be missing)

    crystal_list = entry_data.get("exptl_crystal_grow", [])
    crystal = crystal_list[0] if crystal_list else {}

    print(f"PDB ID: {pdb_id}")
    print("  Crystallization experiment details:")
    if crystal:
        for key, value in crystal.items():
            print(f"    - {key}: {value}")
    else:
        print("    - No crystallization details reported")

    # extract actual PDB ID
    # ---- PubMed / DOI extraction ----

    def extract_pubmed_and_doi(entry_data):
        pubmed = entry_data.get("rcsb_primary_citation", {}).get("pdbx_database_id_pub_med")
        doi = entry_data.get("rcsb_primary_citation", {}).get("pdbx_database_id_doi")

        return pubmed, doi

    pubmed, doi= extract_pubmed_and_doi(entry_data)

    print(f"  PubMed IDs: {pubmed if pubmed else 'Not available'}")
    print(f"  DOIs: {doi if doi else 'Not available'}")
    print("-" * 40)
