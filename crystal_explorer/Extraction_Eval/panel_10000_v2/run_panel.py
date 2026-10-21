"""Run Crystal Explorer's real pipeline (stage 1: RCSB search + metadata,
stage 2: compound extraction) on every target in targets.json, with the
same thresholds as the earlier panels (identity 0.5, E-value 1.0, 30 hits).

    python run_panel.py --shard 0 --nshards 3   # one of 3 parallel workers

Each shard handles every nshards-th target and writes its own
results_shard<k>.json, so shards never write the same file. Resumable: a
target already marked "ok" (in any shard file) is skipped. A target that
takes longer than TIMEOUT_S -- typically one whose homologs are huge
assemblies such as whole ribosomes, with very large mmCIF files -- is
abandoned and logged as a timeout instead of stalling the panel."""
import argparse
import glob
import json
import os
import signal
import sys
import time

import pandas as pd

HERE = os.path.dirname(os.path.abspath(__file__))
PIPELINE_DIR = os.path.join(os.path.dirname(os.path.dirname(HERE)), "pipeline")
sys.path.insert(0, PIPELINE_DIR)

import compound_extraction as ce  # noqa: E402
import pdb_sequence_search as pss  # noqa: E402

IDENTITY = 0.5
EVALUE = 1.0
MAX_HITS = 30
TIMEOUT_S = 300

RUNS_DIR = os.path.join(HERE, "runs")


class TargetTimeout(Exception):
    pass


def _on_alarm(signum, frame):
    raise TargetTimeout(f"timed out after {TIMEOUT_S}s")


def load_all_results() -> dict:
    """Every finished target across the legacy results.json and all shard files."""
    results = {}
    for path in [os.path.join(HERE, "results.json")] + sorted(glob.glob(os.path.join(HERE, "results_shard*.json"))):
        if os.path.exists(path):
            with open(path) as f:
                for r in json.load(f):
                    if r["status"] == "ok" or r["name"] not in results:
                        results[r["name"]] = r
    return results


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--shard", type=int, default=0)
    ap.add_argument("--nshards", type=int, default=1)
    args = ap.parse_args()

    with open(os.path.join(HERE, "targets.json")) as f:
        targets = list(json.load(f).items())
    done = {n for n, r in load_all_results().items() if r["status"] == "ok"}
    mine = [(i, n, info) for i, (n, info) in enumerate(targets, 1)
            if i % args.nshards == args.shard and n not in done]

    shard_path = os.path.join(HERE, f"results_shard{args.shard}.json")
    shard_results = {}
    if os.path.exists(shard_path):
        with open(shard_path) as f:
            shard_results = {r["name"]: r for r in json.load(f)}

    signal.signal(signal.SIGALRM, _on_alarm)
    t_start = time.time()
    for k, (i, name, info) in enumerate(mine, 1):
        protein_dir = os.path.join(RUNS_DIR, name)
        compounds_csv = os.path.join(protein_dir, "Output_compounds.csv")
        os.makedirs(protein_dir, exist_ok=True)
        print(f"=== shard {args.shard}: [{k}/{len(mine)}] target #{i} {name} "
              f"-- {time.time() - t_start:.0f}s elapsed", flush=True)
        t0 = time.time()
        signal.alarm(TIMEOUT_S)
        try:
            rows = pss.find_homologs_with_conditions(
                info["sequence"], identity_cutoff=IDENTITY, evalue_cutoff=EVALUE,
                sequence_type="protein", max_hits=MAX_HITS, sleep=0.1,
                max_workers=8, verbose=False,
            )
            output_csv = os.path.join(protein_dir, "Output.csv")
            pss.filter_experimental_conditions(rows, output_csv=output_csv)
            ce.process_csv(output_csv, compounds_csv)
            n = len(pd.read_csv(output_csv)) if os.path.exists(output_csv) else 0
            status = "ok"
        except Exception as e:  # one bad or oversized target must not stop the panel
            n, status = 0, f"error: {e}"
        finally:
            signal.alarm(0)
        shard_results[name] = {"name": name, "entity": info["entity"], "retained_hits": n,
                               "runtime_s": round(time.time() - t0, 2), "status": status}
        print(shard_results[name], flush=True)
        with open(shard_path, "w") as f:
            json.dump(list(shard_results.values()), f, indent=2)

    ok = [r for r in shard_results.values() if r["status"] == "ok"]
    print(f"\nshard {args.shard} DONE in {time.time() - t_start:.0f}s: {len(ok)} targets ok, "
          f"{sum(r['retained_hits'] for r in ok)} retained hits")


if __name__ == "__main__":
    main()
