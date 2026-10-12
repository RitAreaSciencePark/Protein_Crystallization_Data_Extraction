"""Draw the final, genuinely disjoint sample from the 100-protein held-out
panel's 2,265-hit pool: everything NOT already used in the 100-entry
(60+40), 500-entry, or 1000-entry draws. Only 665 entries remain unused
after those three panels, so this draws all of them rather than a
sub-sample -- it exhausts the pool, meaning every retained hit from the
100-protein panel will have been through manual ground-truth annotation
across the four panels combined.

Unlike the earlier build_sample_*.py scripts (which carried over each
entry's ORIGINAL, now-stale "extracted_compound_string" from the pipeline
run that first produced Output_compounds.csv), this one re-runs the
CURRENT compound_extraction.py against each entry's raw pdbx_details text,
so the attached extraction reflects every dictionary/parsing fix made
since -- a more useful, accurate starting point for the annotator."""
import os
import sys
import json

import pandas as pd

HERE = os.path.dirname(os.path.abspath(__file__))
REPO_ROOT = os.path.dirname(os.path.dirname(os.path.dirname(HERE)))
PIPELINE_DIR = os.path.join(REPO_ROOT, "crystal_explorer", "pipeline")
sys.path.insert(0, PIPELINE_DIR)
import compound_extraction as ce

EVAL_DIR = os.path.join(HERE, "runs")

with open(os.path.join(HERE, "heldout_results.json")) as f:
    results = json.load(f)

already_used = set()
for fname in (
    "heldout_sample_raw.json",
    "heldout_sample_extra40_raw.json",
    "heldout_sample_500_raw.json",
    "heldout_sample_1000_raw.json",
):
    path = os.path.join(HERE, fname)
    if not os.path.exists(path):
        continue
    with open(path) as f:
        prev = json.load(f)
    for e in prev:
        already_used.add((e["protein"], str(e["pdb_id"]), e["pdbx_details"]))

sample = []
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
        clauses = ce.split_into_clauses(str(details))
        matches = []
        for c in clauses:
            matches.extend(ce.match_clause_all(c))
        sample.append({
            "protein": r["name"],
            "pdb_id": pdb_id,
            "pdbx_details": str(details),
            "extracted_compound_string": ce.format_compounds(matches),
        })

print("Final disjoint sample size (should be 665):", len(sample))
with open(os.path.join(HERE, "heldout_sample_665_raw.json"), "w") as f:
    json.dump(sample, f, indent=2)
print("Wrote heldout_sample_665_raw.json with current-code extraction attached.")
