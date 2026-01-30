#!/usr/bin/env python3

import sys
import re
import pandas as pd
import requests

# -------------------------
# Simple reagent vocabulary
# -------------------------

KNOWN_REAGENTS = [
    "glycerol", "hepes", "tris", "mes",
    "ethylene glycol", "acetate",
    "sodium chloride", "magnesium chloride",
    "ammonium sulfate"
]

NAME_TO_CCD = {
    "glycerol": ["GOL"],
    "hepes": ["HEP"],
    "tris": ["TRS"],
    "mes": ["MES"],
    "ethylene glycol": ["EDO"],
    "sodium chloride": ["NA", "CL"],
    "magnesium chloride": ["MG", "CL"],
    "acetate": ["ACT"],
    "ammonium sulfate": ["NH4", "SO4"]
}

ION_CCDS = {"NA", "K", "MG", "CA", "CL", "ZN"}

# -------------------------
# PDBeChem access
# -------------------------

def fetch_pdbechem(ccd):
    url = f"https://www.ebi.ac.uk/pdbe/api/pdb/compound/summary/{ccd}"
    r = requests.get(url, timeout=10)
    if not r.ok:
        return {}
    data = r.json()
    return data.get(ccd.lower(), [{}])[0]

# -------------------------
# Text processing
# -------------------------

def extract_reagents(text):
    text = str(text).lower()
    hits = set()

    for r in KNOWN_REAGENTS:
        if r in text:
            hits.add(r)

    # PEG variants: PEG 400, PEG3350, etc.
    peg_hits = re.findall(r"peg\s*\d*", text)
    hits.update(peg_hits)

    return list(hits)

def reagent_to_ccds(reagent):
    if reagent.startswith("peg"):
        return ["PEG"]
    return NAME_TO_CCD.get(reagent, [])

# -------------------------
# Classification
# -------------------------

def classify_reagent(ccd, chem):
    name = chem.get("name", "").lower()
    mw = chem.get("formula_weight", 0) or 0

    if ccd in ION_CCDS:
        return "ion"
    if ccd in {"HEP", "TRS", "MES"}:
        return "buffer"
    if ccd in {"GOL", "EDO"}:
        return "cryoprotectant"
    if ccd == "PEG" or mw > 500:
        return "precipitant"

    return "additive"

# -------------------------
# Main
# -------------------------

def main(csv_file):
    df = pd.read_csv(csv_file)

    if not {"pdb_id", "pdbx_details"}.issubset(df.columns):
        sys.exit("ERROR: CSV must contain pdb_id and pdbx_details columns")

    rows = []

    for _, row in df.iterrows():
        pdb_id = row["pdb_id"]
        details = row["pdbx_details"]

        reagents = extract_reagents(details)

        for reagent in reagents:
            ccds = reagent_to_ccds(reagent)

            for ccd in ccds:
                chem = fetch_pdbechem(ccd)
                rows.append({
                    "pdb_id": pdb_id,
                    "reagent_text": reagent,
                    "ccd": ccd,
                    "chemical_name": chem.get("name"),
                    "formula_weight": chem.get("formula_weight"),
                    "class": classify_reagent(ccd, chem)
                })

    out = pd.DataFrame(rows)
    out.to_csv(sys.stdout, index=False)

if __name__ == "__main__":
    if len(sys.argv) != 2:
        sys.exit(f"Usage: {sys.argv[0]} input.csv")
    main(sys.argv[1])
