"""Fetch the structural context of every entry of the held-out panel, for the
context-controlled version of the identity-vs-conditions analysis
(context_analysis.py): the oligomeric state of the first biological
assembly, the chain count of each polymer entity, and the non-polymer
components (ligands, ions, additives).

Sequences and UniProt accessions per polymer entity are already in
entry_metadata.json (fetch_metadata.py); this adds only what that file lacks.

    python fetch_context.py      # -> entry_context.json (resumable)
"""
import json
import os
import time

import requests

HERE = os.path.dirname(os.path.abspath(__file__))
META = os.path.join(HERE, "entry_metadata.json")
OUT = os.path.join(HERE, "entry_context.json")
URL = "https://data.rcsb.org/graphql"
BATCH = 200

QUERY = """query($ids:[String!]!){entries(entry_ids:$ids){
  rcsb_id
  polymer_entities{rcsb_id rcsb_polymer_entity_container_identifiers{auth_asym_ids}}
  assemblies{rcsb_id rcsb_struct_symmetry{kind oligomeric_state stoichiometry}}
  nonpolymer_entities{nonpolymer_comp{chem_comp{id}}}
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
    chains = {p["rcsb_id"]: len((p.get("rcsb_polymer_entity_container_identifiers") or {})
                                .get("auth_asym_ids") or [])
              for p in e.get("polymer_entities") or []}
    asm = sorted(e.get("assemblies") or [], key=lambda a: a["rcsb_id"])
    olig, stoich = None, None
    if asm:
        for s in asm[0].get("rcsb_struct_symmetry") or []:
            if s.get("kind") == "Global Symmetry":
                olig, stoich = s.get("oligomeric_state"), s.get("stoichiometry")
                break
    ligands = sorted({n["nonpolymer_comp"]["chem_comp"]["id"]
                      for n in e.get("nonpolymer_entities") or []
                      if (n.get("nonpolymer_comp") or {}).get("chem_comp")})
    return {"chains_per_entity": chains, "assembly1_oligomeric_state": olig,
            "assembly1_stoichiometry": stoich, "nonpolymer": ligands}


def main():
    with open(META) as f:
        ids = sorted(json.load(f))
    ctx = {}
    if os.path.exists(OUT):
        with open(OUT) as f:
            ctx = json.load(f)
    todo = [i for i in ids if i not in ctx]
    print(f"{len(ids)} entries, {len(todo)} to fetch", flush=True)
    for k in range(0, len(todo), BATCH):
        chunk = todo[k:k + BATCH]
        for e in fetch(chunk):
            ctx[e["rcsb_id"]] = simplify(e)
        with open(OUT, "w") as f:
            json.dump(ctx, f)
        print(f"  {min(k + BATCH, len(todo))}/{len(todo)}", flush=True)
    missing = [i for i in ids if i not in ctx]
    print(f"done: {len(ctx)} fetched, {len(missing)} missing {missing[:10]}")


if __name__ == "__main__":
    main()
