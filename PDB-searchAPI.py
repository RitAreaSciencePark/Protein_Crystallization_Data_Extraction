import requests
import json

ROW_COUNT = 10

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
data = requests.get(f'https://data.rcsb.org/rest/v1/core/entry/{identifier}').json()

#Print all top-level data fields (preview)
print("We get a lot of data (and still, is not everything): ")
for labels in data:
  print("\t - "+ labels + ": " + str(data[labels])[:50] + ".......")
print('\n Let\'s focus on the last bit of information:')

#print all external references (e.g., doi, pubmed, etc.)
#print(data['rcsb_external_references'])

# Loop over all results and print some information about the experimental conditions
search_data = result.json()

if search_data["total_count"] == 0:
    raise ValueError("No matching PDB entries found")

print("\nExtracting experimental conditions:\n")

for hit in search_data["result_set"]:
    pdb_id = hit["identifier"]

    # Fetch entry metadata
    entry_data = requests.get(
        f"https://data.rcsb.org/rest/v1/core/entry/{pdb_id}"
    ).json()

    # Experimental method
    method = entry_data.get("exptl", [{}])[0].get("method", "NA")

    # Crystallization conditions (may be missing)
    crystal_info = entry_data.get("exptl_crystal_grow", [{}])
    crystal_info = crystal_info[0] if crystal_info else {}

    ph = crystal_info.get("pH", "NA")
    temperature = crystal_info.get("temp", "NA")

    print(f"PDB ID: {pdb_id}")
    print(f"  Method: {method}")
    print(f"  pH: {ph}")
    print(f"  Temperature: {temperature}")
    print("-" * 40)
