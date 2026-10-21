"""Build the 10,000-entry panel and its 1,000-entry audit sample.

  * Pool: every retained hit across runs/<target>/Output.csv with a
    non-empty pdbx_details, one row per PDB ID (the first target that
    retrieved it wins), excluding every PDB ID seen in earlier panels.
  * panel_10000_raw.json: 10,000 entries drawn at random (seed 2026) from
    that pool, each with the CURRENT extractor's output attached.
  * audit_1000_raw.json + batches/batch_NN_raw.json: a random 1,000-entry
    subset (seed 2027) of the panel, split into 20 batches of 50 for the
    manual ground-truth audit.

The sha256 of compound_extraction.py is stored in panel_meta.json: the
audit scores that exact version, so the extractor must not be changed
while the audit is in progress."""
import hashlib
import json
import os
import random
import sys

import pandas as pd

HERE = os.path.dirname(os.path.abspath(__file__))
PIPELINE_DIR = os.path.join(os.path.dirname(os.path.dirname(HERE)), "pipeline")
sys.path.insert(0, PIPELINE_DIR)
sys.path.insert(0, HERE)

import compound_extraction as ce  # noqa: E402
from run_panel import load_all_results  # noqa: E402
from select_targets import previously_used_pdb_ids  # noqa: E402

PANEL_SIZE = 10_000
AUDIT_SIZE = 1_000
BATCH_SIZE = 50


def extract(details: str) -> str:
    return ce.format_compounds([m for c in ce.split_into_clauses(details) for m in ce.match_clause_all(c)])


def main():
    results = [r for r in load_all_results().values() if r["status"] == "ok"]
    used = previously_used_pdb_ids()

    pool, seen = [], set()
    for r in results:
        csv = os.path.join(HERE, "runs", r["name"], "Output.csv")
        if not os.path.exists(csv):
            continue
        df = pd.read_csv(csv)
        id_col = next(c for c in df.columns if c.lower() == "pdb_id")
        for _, row in df.iterrows():
            pdb_id = str(row[id_col]).upper()
            details = row.get("pdbx_details")
            if pdb_id in seen or pdb_id in used or pd.isna(details) or not str(details).strip():
                continue
            seen.add(pdb_id)
            pool.append({"protein": r["name"], "pdb_id": pdb_id, "pdbx_details": str(details)})
    print(f"pool: {len(pool)} unique, unused PDB entries with pdbx_details")
    if len(pool) < PANEL_SIZE:
        sys.exit(f"pool too small for a {PANEL_SIZE}-entry panel -- run more targets")

    panel = random.Random(2028).sample(pool, PANEL_SIZE)
    for e in panel:
        e["extracted_compound_string"] = extract(e["pdbx_details"])
    audit = random.Random(2029).sample(panel, AUDIT_SIZE)

    with open(os.path.join(HERE, "panel_10000_raw.json"), "w") as f:
        json.dump(panel, f, indent=2)
    with open(os.path.join(HERE, "audit_1000_raw.json"), "w") as f:
        json.dump(audit, f, indent=2)
    os.makedirs(os.path.join(HERE, "batches"), exist_ok=True)
    for i in range(0, AUDIT_SIZE, BATCH_SIZE):
        with open(os.path.join(HERE, "batches", f"batch_{i // BATCH_SIZE + 1:02d}_raw.json"), "w") as f:
            json.dump(audit[i:i + BATCH_SIZE], f, indent=2)

    with open(os.path.join(PIPELINE_DIR, "compound_extraction.py"), "rb") as f:
        sha = hashlib.sha256(f.read()).hexdigest()
    with open(os.path.join(HERE, "panel_meta.json"), "w") as f:
        json.dump({"pool_size": len(pool), "panel_size": PANEL_SIZE, "audit_size": AUDIT_SIZE,
                   "panel_seed": 2028, "audit_seed": 2029, "targets_ok": len(results),
                   "compound_extraction_sha256": sha}, f, indent=2)
    print(f"wrote panel ({PANEL_SIZE}), audit sample ({AUDIT_SIZE}), "
          f"{AUDIT_SIZE // BATCH_SIZE} batches; extractor sha256 {sha[:12]}")


if __name__ == "__main__":
    main()
