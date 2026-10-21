"""Select the target proteins for the SECOND 10,000-entry evaluation panel
(after the fixes made from panel_10000). Same method as panel_10000, but
seed 2028, and every PDB ID from panel_10000's runs and every cluster it
selected are excluded too.

Original description:
Select the target proteins for the 10,000-entry evaluation panel.

Instead of hand-typing ~800 (name, PDB ID) pairs -- the 100-protein panel
already hit stale/wrong remembered IDs -- targets are drawn reproducibly
from RCSB itself: X-ray protein polymer entities (60-900 aa) are grouped
into 30%-sequence-identity clusters, and one representative is taken per
randomly sampled cluster (seed 2026). Clusters at 30% identity are mutually
unrelated, so the targets are diverse by construction.

Two filters keep the panel genuinely new:
  * only clusters with >= MIN_CLUSTER_SIZE members, so a 30-hit homology
    search has enough structures to return;
  * any cluster containing a PDB ID already seen in an earlier panel's
    pipeline runs (../runs, ../heldout/runs) is dropped, so no protein
    family from the development or earlier held-out panels comes back.

Writes targets.json (name, representative PDB/entity, sequence) and
selection_log.json next to this script."""
import glob
import json
import os
import random
import time

import pandas as pd
import requests

HERE = os.path.dirname(os.path.abspath(__file__))
EVAL_DIR = os.path.dirname(HERE)
SEARCH_URL = "https://search.rcsb.org/rcsbsearch/v2/query"
SEED = 2028
N_TARGETS = 800
MIN_CLUSTER_SIZE = 20
PAGE = 1000


def previously_used_pdb_ids() -> set:
    ids = set()
    for csv in glob.glob(os.path.join(EVAL_DIR, "runs", "*", "Output.csv")) + \
               glob.glob(os.path.join(EVAL_DIR, "heldout", "runs", "*", "Output.csv")) + \
               glob.glob(os.path.join(EVAL_DIR, "panel_10000", "runs", "*", "Output.csv")):
        df = pd.read_csv(csv, usecols=lambda c: c.lower() == "pdb_id")
        ids.update(str(x).upper() for x in df.iloc[:, 0].dropna())
    return ids


def fetch_clusters() -> list:
    query = {
        "query": {"type": "group", "logical_operator": "and", "nodes": [
            {"type": "terminal", "service": "text", "parameters": {
                "attribute": "exptl.method", "operator": "exact_match", "value": "X-RAY DIFFRACTION"}},
            {"type": "terminal", "service": "text", "parameters": {
                "attribute": "entity_poly.rcsb_entity_polymer_type", "operator": "exact_match", "value": "Protein"}},
            {"type": "terminal", "service": "text", "parameters": {
                "attribute": "entity_poly.rcsb_sample_sequence_length", "operator": "range",
                "value": {"from": 60, "to": 900, "include_lower": True, "include_upper": True}}},
        ]},
        "return_type": "polymer_entity",
        "request_options": {
            "group_by": {"aggregation_method": "sequence_identity", "similarity_cutoff": 30},
            "group_by_return_type": "groups",
            "paginate": {"start": 0, "rows": PAGE},
        },
    }
    clusters, start = [], 0
    while True:
        query["request_options"]["paginate"]["start"] = start
        r = requests.post(SEARCH_URL, json=query, timeout=120)
        r.raise_for_status()
        groups = r.json().get("group_set", [])
        if not groups:
            break
        for g in groups:
            clusters.append({
                "cluster": g["identifier"],
                "count": g["count"],
                "representative": g["representative"]["identifier"],
                "members": [m["identifier"].split("_")[0].upper() for m in g["result_set"]],
            })
        print(f"fetched {len(clusters)} clusters", flush=True)
        start += PAGE
        time.sleep(0.2)
    return clusters


def fetch_entity(entity_id: str):
    entry, ent = entity_id.split("_")
    r = requests.get(f"https://data.rcsb.org/rest/v1/core/polymer_entity/{entry}/{ent}", timeout=30)
    r.raise_for_status()
    d = r.json()
    seq = d["entity_poly"]["pdbx_seq_one_letter_code_can"].replace("\n", "")
    desc = (d.get("rcsb_polymer_entity") or {}).get("pdbx_description") or entity_id
    return seq, desc


def safe_name(desc: str, entity_id: str) -> str:
    base = "".join(ch if ch.isalnum() else "_" for ch in desc).strip("_")[:50] or "protein"
    return f"{base}__{entity_id}"


def main():
    used = previously_used_pdb_ids()
    print(f"{len(used)} PDB IDs already used in earlier panels")
    clusters = fetch_clusters()
    with open(os.path.join(EVAL_DIR, "panel_10000", "targets.json")) as f:
        used_clusters = {t["cluster"] for t in json.load(f).values()}
    eligible = [c for c in clusters
                if c["count"] >= MIN_CLUSTER_SIZE and not used.intersection(c["members"])
                and c["cluster"] not in used_clusters]
    print(f"{len(clusters)} clusters, {len(eligible)} eligible (size >= {MIN_CLUSTER_SIZE}, no overlap)")

    random.Random(SEED).shuffle(eligible)
    targets, failed = {}, []
    for c in eligible:
        if len(targets) >= N_TARGETS:
            break
        rep = c["representative"]
        try:
            seq, desc = fetch_entity(rep)
        except Exception as e:  # stale entity, network hiccup: skip and keep going
            failed.append({"entity": rep, "reason": str(e)})
            continue
        if len(seq) < 30:
            failed.append({"entity": rep, "reason": f"sequence too short ({len(seq)})"})
            continue
        name = safe_name(desc, rep)
        targets[name] = {"pdb_id": rep.split("_")[0], "entity": rep, "description": desc,
                         "cluster": c["cluster"], "cluster_size": c["count"],
                         "sequence": seq, "length": len(seq)}
        print(f"[{len(targets)}/{N_TARGETS}] {name} ({len(seq)} aa)", flush=True)
        time.sleep(0.05)

    with open(os.path.join(HERE, "targets.json"), "w") as f:
        json.dump(targets, f, indent=2)
    with open(os.path.join(HERE, "selection_log.json"), "w") as f:
        json.dump({"seed": SEED, "n_clusters": len(clusters), "n_eligible": len(eligible),
                   "min_cluster_size": MIN_CLUSTER_SIZE, "previously_used_ids": len(used),
                   "n_targets": len(targets), "failed": failed}, f, indent=2)
    print(f"\n{len(targets)} targets written, {len(failed)} failed")


if __name__ == "__main__":
    main()
