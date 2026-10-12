"""Draw a fresh, independent 1000-entry sample for manual ground-truth
annotation, from the same 100-protein held-out panel's 2,265-hit pool,
explicitly excluding every entry already used in the earlier 100-entry
(60+40) and 500-entry samples, so this is a genuinely new, disjoint
draw."""
import os
import json
import random

import pandas as pd

HERE = os.path.dirname(os.path.abspath(__file__))
EVAL_DIR = os.path.join(HERE, "runs")
SAMPLE_SIZE = 1000
SEED = 1000

with open(os.path.join(HERE, "heldout_results.json")) as f:
    results = json.load(f)

already_used = set()
for fname in (
    "heldout_sample_raw.json",
    "heldout_sample_extra40_raw.json",
    "heldout_sample_500_raw.json",
):
    path = os.path.join(HERE, fname)
    if not os.path.exists(path):
        continue
    with open(path) as f:
        prev = json.load(f)
    for e in prev:
        already_used.add((e["protein"], str(e["pdb_id"]), e["pdbx_details"]))

pool = []
for r in results:
    if r.get("status") != "ok":
        continue
    csv_path = os.path.join(EVAL_DIR, r["name"], "Output_compounds.csv")
    if not os.path.exists(csv_path):
        continue
    df = pd.read_csv(csv_path)
    for _, row in df.iterrows():
        details = row.get("pdbx_details")
        if pd.isna(details) or not str(details).strip():
            continue
        pdb_id = row.get("pdb_id", row.get("PDB_ID"))
        key = (r["name"], str(pdb_id), str(details))
        if key in already_used:
            continue
        pool.append({
            "protein": r["name"],
            "pdb_id": pdb_id,
            "pdbx_details": str(details),
            "extracted_compound_string": row.get("compound") if pd.notna(row.get("compound")) else "",
        })

print("Pool size (unused retained hits with pdbx_details):", len(pool))
random.seed(SEED)
sample = random.sample(pool, min(SAMPLE_SIZE, len(pool)))
with open(os.path.join(HERE, "heldout_sample_1000_raw.json"), "w") as f:
    json.dump(sample, f, indent=2)
print(f"Wrote {len(sample)}-entry fresh sample for manual ground-truth annotation.")
