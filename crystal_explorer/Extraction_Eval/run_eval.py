"""Run Crystal Explorer's real pipeline (stage 1: search+metadata, stage 2:
compound extraction) on the 8 benchmark proteins from eval_targets.json.
Produces one Output_compounds.csv per protein under runs/, which
build_ground_truth_sample.py then draws a random sample from for manual
annotation. Must be run with a Python that has crystal_explorer's
requirements.txt installed (e.g. crystal_explorer/.venv/bin/python3
Extraction_Eval/run_eval.py)."""
import sys, os, json, time

HERE = os.path.dirname(os.path.abspath(__file__))
REPO_ROOT = os.path.dirname(os.path.dirname(HERE))
PIPELINE_DIR = os.path.join(REPO_ROOT, "crystal_explorer", "pipeline")
sys.path.insert(0, PIPELINE_DIR)

import pdb_sequence_search as pss
import compound_extraction as ce
import pandas as pd

with open(os.path.join(HERE, "eval_targets.json")) as f:
    targets = json.load(f)

EVAL_DIR = os.path.join(HERE, "runs")
os.makedirs(EVAL_DIR, exist_ok=True)

IDENTITY = 0.5
EVALUE = 1.0
MAX_HITS = 30

results = []
for name, info in targets.items():
    seq = info["sequence"]
    protein_dir = os.path.join(EVAL_DIR, name)
    os.makedirs(protein_dir, exist_ok=True)

    print(f"\n=== {name} ({info['pdb_id']}, {info['length']} aa) ===", flush=True)
    t0 = time.time()
    rows = pss.find_homologs_with_conditions(
        seq, identity_cutoff=IDENTITY, evalue_cutoff=EVALUE,
        sequence_type="protein", max_hits=MAX_HITS, sleep=0.1,
        max_workers=8, verbose=False,
    )
    output_csv = os.path.join(protein_dir, "Output.csv")
    pss.filter_experimental_conditions(rows, output_csv=output_csv)
    t1 = time.time()

    compounds_csv = os.path.join(protein_dir, "Output_compounds.csv")
    ce.process_csv(output_csv, compounds_csv)

    df = pd.read_csv(output_csv)
    result = {"name": name, "pdb_id": info["pdb_id"], "retained_hits": len(df),
              "runtime_s": round(time.time() - t0, 2)}
    print(result, flush=True)
    results.append(result)

with open(os.path.join(HERE, "eval_results.json"), "w") as f:
    json.dump(results, f, indent=2)
print("\nDONE")
