import requests
import json

ROW_COUNT = 10

#this query is to navigate the PDB database and find the proteins base on the sequences determined by X-ray diffraction ,methods

query ={
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
            "value": "ELECTRON MICROSCOPY"
          }
        }
      ]
    },
    "request_options": {
      "paginate": {
        "start": 0,
        "rows": ROW_COUNT
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
