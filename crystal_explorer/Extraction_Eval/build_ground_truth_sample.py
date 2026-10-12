"""Pool every retained hit's pdbx_details + compound_extraction.py's own
extraction across all 8 proteins' runs/ output, then draw a random sample
for manual ground-truth annotation. Writes
ground_truth_sample_raw.json (unannotated) next to this script; the
paper's authors then read each entry by hand and produced
ground_truth_annotations.json (tp/fp/fn judgments against the source
text), which analyze_ground_truth.py scores."""
import os
import json
import random

import pandas as pd

HERE = os.path.dirname(os.path.abspath(__file__))
EVAL_DIR = os.path.join(HERE, "runs")
SAMPLE_SIZE = 60
SEED = 7

with open(os.path.join(HERE, "eval_results.json")) as f:
    results = json.load(f)

pool = []
for r in results:
    df = pd.read_csv(os.path.join(EVAL_DIR, r["name"], "Output_compounds.csv"))
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

print("Pool size (retained hits with pdbx_details):", len(pool))
random.seed(SEED)
sample = random.sample(pool, SAMPLE_SIZE)
with open(os.path.join(HERE, "ground_truth_sample_raw.json"), "w") as f:
    json.dump(sample, f, indent=2)
print(f"Wrote {SAMPLE_SIZE}-entry sample for manual ground-truth annotation.")
