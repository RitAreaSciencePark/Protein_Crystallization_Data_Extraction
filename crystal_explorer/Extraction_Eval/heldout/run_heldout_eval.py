"""Run Crystal Explorer's real pipeline (stage 1: search+metadata, stage 2:
compound extraction) on the 100-protein held-out validation panel. Same
code path, same thresholds, as the original 8-protein development panel's
run_eval.py -- this is a fresh, independent set of proteins never seen
during the five rounds of dictionary/parsing fixes, to check whether those
fixes generalize rather than just fitting the development sample."""
import sys, os, json, time

HERE = os.path.dirname(os.path.abspath(__file__))
REPO_ROOT = os.path.dirname(os.path.dirname(os.path.dirname(HERE)))
PIPELINE_DIR = os.path.join(REPO_ROOT, "crystal_explorer", "pipeline")
sys.path.insert(0, PIPELINE_DIR)

import pdb_sequence_search as pss
import compound_extraction as ce
import pandas as pd

with open(os.path.join(HERE, "heldout_targets.json")) as f:
    targets = json.load(f)

EVAL_DIR = os.path.join(HERE, "runs")
os.makedirs(EVAL_DIR, exist_ok=True)

IDENTITY = 0.5
EVALUE = 1.0
MAX_HITS = 30

results = []
results_path = os.path.join(HERE, "heldout_results.json")
t_start = time.time()
for i, (name, info) in enumerate(targets.items(), 1):
    seq = info["sequence"]
    protein_dir = os.path.join(EVAL_DIR, name)
    os.makedirs(protein_dir, exist_ok=True)

    print(f"\n=== [{i}/{len(targets)}] {name} ({info['pdb_id']}, {info['length']} aa) "
          f"-- {time.time()-t_start:.0f}s elapsed ===", flush=True)
    t0 = time.time()
    try:
        rows = pss.find_homologs_with_conditions(
            seq, identity_cutoff=IDENTITY, evalue_cutoff=EVALUE,
            sequence_type="protein", max_hits=MAX_HITS, sleep=0.1,
            max_workers=8, verbose=False,
        )
        output_csv = os.path.join(protein_dir, "Output.csv")
        pss.filter_experimental_conditions(rows, output_csv=output_csv)

        compounds_csv = os.path.join(protein_dir, "Output_compounds.csv")
        ce.process_csv(output_csv, compounds_csv, use_llm_fallback=False)

        df = pd.read_csv(output_csv)
        result = {"name": name, "pdb_id": info["pdb_id"], "retained_hits": len(df),
                  "runtime_s": round(time.time() - t0, 2), "status": "ok"}
    except Exception as e:
        result = {"name": name, "pdb_id": info["pdb_id"], "retained_hits": 0,
                   "runtime_s": round(time.time() - t0, 2), "status": f"error: {e}"}
    print(result, flush=True)
    results.append(result)

    # Write incrementally so partial progress survives if this is interrupted.
    with open(results_path, "w") as f:
        json.dump(results, f, indent=2)

print(f"\nDONE in {time.time()-t_start:.0f}s. Total retained hits: "
      f"{sum(r['retained_hits'] for r in results)}")
