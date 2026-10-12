"""Fetch the 8 benchmark query sequences (deposited chain, longest polymer
record) from RCSB for the ACDSA 2027 paper's Evaluation section. Run from
anywhere; writes eval_targets.json next to this script."""
import os
import json
import requests

HERE = os.path.dirname(os.path.abspath(__file__))

targets = {
    "lysozyme_hen":   "6LYZ",
    "lysozyme_T4":    "2LZM",
    "RNaseA_bovine":  "7RSA",
    "trypsin_bovine": "3PTN",
    "Mpro_SARS2":     "6LU7",
    "GFP":            "1EMA",
    "myoglobin_sperm_whale": "1MBN",
    "carbonic_anhydrase_II_human": "1CA2",
}

out = {}
for name, pdb_id in targets.items():
    url = f"https://www.rcsb.org/fasta/entry/{pdb_id}"
    r = requests.get(url, timeout=20)
    r.raise_for_status()
    text = r.text.strip()
    records = [rec for rec in text.split(">") if rec.strip()]
    parsed = []
    for rec in records:
        lines = rec.strip().split("\n")
        header = lines[0]
        seq = "".join(lines[1:]).strip()
        parsed.append((header, seq))
    header, seq = max(parsed, key=lambda hs: len(hs[1]))
    out[name] = {"pdb_id": pdb_id, "header": header, "sequence": seq, "length": len(seq)}
    print(name, pdb_id, len(seq), header[:80])

with open(os.path.join(HERE, "eval_targets.json"), "w") as f:
    json.dump(out, f, indent=2)
