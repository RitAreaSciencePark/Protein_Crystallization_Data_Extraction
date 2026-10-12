"""Pool every retained hit's pdbx_details + compound_extraction.py's own
extraction across all 100 held-out proteins' runs/ output, then draw a
random sample for manual ground-truth annotation -- same methodology as
the original 8-protein panel's build_ground_truth_sample.py, applied to
this independent, never-seen-during-development panel."""
import os
import json
import random

import pandas as pd

HERE = os.path.dirname(os.path.abspath(__file__))
EVAL_DIR = os.path.join(HERE, "runs")
SAMPLE_SIZE = 60
SEED = 11

with open(os.path.join(HERE, "heldout_results.json")) as f:
    results = json.load(f)

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
        pool.append({
            "protein": r["name"],
            "pdb_id": row.get("pdb_id", row.get("PDB_ID")),
            "pdbx_details": str(details),
            "extracted_compound_string": row.get("compound") if pd.notna(row.get("compound")) else "",
        })

print("Pool size (retained hits with pdbx_details, held-out panel):", len(pool))
random.seed(SEED)
sample = random.sample(pool, min(SAMPLE_SIZE, len(pool)))
with open(os.path.join(HERE, "heldout_sample_raw.json"), "w") as f:
    json.dump(sample, f, indent=2)
print(f"Wrote {len(sample)}-entry sample for manual ground-truth annotation.")
