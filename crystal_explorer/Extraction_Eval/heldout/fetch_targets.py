"""Fetch real deposited sequences for the held-out validation panel from
RCSB's FASTA endpoint (same method as the original 8-protein panel's
fetch_seqs.py). Keeps the first 100 candidates that resolve to a real
PDB entry with a usable polymer sequence; skips and logs any that 404 or
otherwise fail (stale/incorrect remembered PDB IDs). Writes
heldout_targets.json next to this script."""
import os
import sys
import json
import time
import requests

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from candidates import CANDIDATES

HERE = os.path.dirname(os.path.abspath(__file__))
TARGET_COUNT = 100

out = {}
failed = []
for name, pdb_id in CANDIDATES:
    if len(out) >= TARGET_COUNT:
        break
    url = f"https://www.rcsb.org/fasta/entry/{pdb_id}"
    try:
        r = requests.get(url, timeout=15)
        if r.status_code != 200:
            failed.append((name, pdb_id, f"HTTP {r.status_code}"))
            continue
        text = r.text.strip()
        records = [rec for rec in text.split(">") if rec.strip()]
        parsed = []
        for rec in records:
            lines = rec.strip().split("\n")
            header = lines[0]
            seq = "".join(lines[1:]).strip()
            parsed.append((header, seq))
        if not parsed:
            failed.append((name, pdb_id, "no FASTA records"))
            continue
        header, seq = max(parsed, key=lambda hs: len(hs[1]))
        if len(seq) < 30:
            failed.append((name, pdb_id, f"sequence too short ({len(seq)} aa)"))
            continue
        out[name] = {"pdb_id": pdb_id, "header": header, "sequence": seq, "length": len(seq)}
        print(f"[{len(out)}/{TARGET_COUNT}] {name} {pdb_id} {len(seq)} aa")
    except Exception as e:
        failed.append((name, pdb_id, str(e)))
    time.sleep(0.1)

print(f"\nResolved {len(out)} targets, {len(failed)} failed.")
if failed:
    print("Failed:")
    for name, pdb_id, reason in failed:
        print(f"  {name} ({pdb_id}): {reason}")

with open(os.path.join(HERE, "heldout_targets.json"), "w") as f:
    json.dump(out, f, indent=2)
with open(os.path.join(HERE, "heldout_failed.json"), "w") as f:
    json.dump(failed, f, indent=2)
