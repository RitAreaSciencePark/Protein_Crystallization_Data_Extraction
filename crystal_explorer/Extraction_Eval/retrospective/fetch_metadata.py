"""Fetch the per-entry metadata the retrospective analysis needs for every
entry of the held-out 10,000-entry panel (../panel_10000_v2/panel_10000_raw.json):
polymer entity sequences and UniProt accessions, crystal-growth pH and
temperature, primary-citation PubMed ID and deposition date.

The panel itself only stores the target group, the PDB ID, pdbx_details and
the extractor output; the per-target pipeline CSVs (which had identity and
UniProt) were deleted to save space, so this pulls the same fields again
from the RCSB Data API in batches.

    python fetch_metadata.py      # -> entry_metadata.json (resumable)
"""
import json
import os
import time

import requests

HERE = os.path.dirname(os.path.abspath(__file__))
PANEL = os.path.join(HERE, "..", "panel_10000_v2", "panel_10000_raw.json")
OUT = os.path.join(HERE, "entry_metadata.json")
URL = "https://data.rcsb.org/graphql"
BATCH = 200

QUERY = """query($ids:[String!]!){entries(entry_ids:$ids){
  rcsb_id
  exptl_crystal_grow{pH temp}
  rcsb_primary_citation{pdbx_database_id_PubMed}
  rcsb_accession_info{deposit_date}
  polymer_entities{
    rcsb_id
    entity_poly{type pdbx_seq_one_letter_code_can}
    rcsb_polymer_entity_container_identifiers{uniprot_ids}
  }
}}"""


def fetch(ids):
    for attempt in range(4):
        try:
            r = requests.post(URL, json={"query": QUERY, "variables": {"ids": ids}}, timeout=120)
            r.raise_for_status()
            return r.json()["data"]["entries"]
        except Exception as e:  # transient RCSB errors
            print(f"  retry {attempt + 1}: {e}", flush=True)
            time.sleep(5 * (attempt + 1))
    raise RuntimeError("RCSB GraphQL failed repeatedly")


def simplify(e):
    grow = e.get("exptl_crystal_grow") or [{}]
    cit = e.get("rcsb_primary_citation") or {}
    ents = []
    for p in e.get("polymer_entities") or []:
        poly = p.get("entity_poly") or {}
        ids = (p.get("rcsb_polymer_entity_container_identifiers") or {}).get("uniprot_ids") or []
        ents.append({"entity": p["rcsb_id"], "type": poly.get("type"),
                     "seq": poly.get("pdbx_seq_one_letter_code_can") or "", "uniprot": ids})
    return {"pH": grow[0].get("pH"), "temp": grow[0].get("temp"),
            "pubmed": cit.get("pdbx_database_id_PubMed"),
            "deposit_date": (e.get("rcsb_accession_info") or {}).get("deposit_date"),
            "entities": ents}


def main():
    with open(PANEL) as f:
        ids = sorted({e["pdb_id"].upper() for e in json.load(f)})
    meta = {}
    if os.path.exists(OUT):
        with open(OUT) as f:
            meta = json.load(f)
    todo = [i for i in ids if i not in meta]
    print(f"{len(ids)} entries, {len(todo)} to fetch", flush=True)
    for k in range(0, len(todo), BATCH):
        chunk = todo[k:k + BATCH]
        for e in fetch(chunk):
            meta[e["rcsb_id"]] = simplify(e)
        with open(OUT, "w") as f:
            json.dump(meta, f)
        print(f"  {min(k + BATCH, len(todo))}/{len(todo)}", flush=True)
    missing = [i for i in ids if i not in meta]
    print(f"done: {len(meta)} fetched, {len(missing)} missing {missing[:10]}")


if __name__ == "__main__":
    main()
