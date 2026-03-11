# extract_structures.py
import sys
import types
import pickle
import pandas as pd
from pathlib import Path
import ast
import re

def extract_structures_with_metadata(structures_file, filtered_csv_file, output_csv_file):
    """
    Extract PDB_ID, COMPOUNDS from pickle and merge SCORE, PH, PUBMED_ID from filtered CSV.
    Also extracts PEG and PEG_CONS from COMPOUNDS into separate columns.
    Saves output CSV.
    """
    # -------------------------------
    # Dummy module/class to allow unpickling
    # -------------------------------
    fake_module = types.ModuleType("pdb_crystal_database")
    sys.modules["pdb_crystal_database"] = fake_module

    class Structure:
        def __init__(self, *args, **kwargs):
            for k, v in kwargs.items():
                setattr(self, k, v)

    fake_module.Structure = Structure

    # -------------------------------
    # Load pickle
    # -------------------------------
    with open(structures_file, "rb") as f:
        structures_list = pickle.load(f)
    print(f"Loaded {len(structures_list)} structures from pickle.")

    if not structures_list:
        print("Warning: Pickle file is empty. Exiting.")
        return

    # -------------------------------
    # Load filtered CSV
    # -------------------------------
    filtered_df = pd.read_csv(filtered_csv_file, on_bad_lines="skip")
    filtered_df["PDB_ID"] = filtered_df["PDB_ID"].astype(str).str.strip().str.upper()

    # -------------------------------
    # Identify relevant columns in filtered CSV
    # -------------------------------
    def find_column(df, keywords):
        for col in df.columns:
            if any(kw.lower() in col.lower() for kw in keywords):
                return col
        return None

    score_col = find_column(filtered_df, ["score"])
    pubmed_col = find_column(filtered_df, ["pubmed"])
    ph_col = find_column(filtered_df, ["ph"])
    ligand_col = find_column(filtered_df, ["ligand"])

    score_dict = dict(zip(filtered_df["PDB_ID"], filtered_df[score_col])) if score_col else {}
    pubmed_dict = dict(zip(filtered_df["PDB_ID"], filtered_df[pubmed_col])) if pubmed_col else {}
    ph_dict = dict(zip(filtered_df["PDB_ID"], filtered_df[ph_col])) if ph_col else {}
    ligand_dict = dict(zip(filtered_df["PDB_ID"], filtered_df[ligand_col])) if ligand_col else {}

    # -------------------------------
    # Detect PDB_ID and COMPOUNDS in pickle
    # -------------------------------
    first_obj = structures_list[0]
    keys = first_obj.__dict__.keys()
    pdb_attr = next((k for k in keys if 'pdb' in k.lower()), None)
    comp_attr = next((k for k in keys if 'compound' in k.lower()), None)

    if not pdb_attr or not comp_attr:
        raise ValueError("Pickle objects must contain PDB_ID and COMPOUNDS attributes.")

    # -------------------------------
    # Helper: extract PEG
    # -------------------------------
    def extract_peg(compounds):
        """
        compounds: list or string representation of list
        Returns: cleaned_compounds, peg, peg_cons
        """
        peg = None
        peg_cons = None
        cleaned_compounds = []

        if not compounds:
            return cleaned_compounds, peg, peg_cons

        # If string, convert to list
        if isinstance(compounds, str):
            try:
                compounds_list = ast.literal_eval(compounds)
            except Exception:
                compounds_list = []
        elif isinstance(compounds, list):
            compounds_list = compounds
        else:
            # Unexpected type
            compounds_list = []

        i = 0
        while i < len(compounds_list):
            compound = compounds_list[i]
            conc = compounds_list[i+1] if i+1 < len(compounds_list) else None

            if compound and re.search(r"\bPEG\b", str(compound), re.IGNORECASE):
                peg = compound
                peg_cons = conc
            else:
                cleaned_compounds.append(compound)
                cleaned_compounds.append(conc)
            i += 2

        return cleaned_compounds, peg, peg_cons

    # -------------------------------
    # Extract and merge data
    # -------------------------------
    data = []
    found_pdb_ids = set()

    for obj in structures_list:
        pdb_val = str(getattr(obj, pdb_attr, "")).strip().upper()
        if pdb_val in filtered_df["PDB_ID"].values:
            comp_val = getattr(obj, comp_attr, "")
            score_val = score_dict.get(pdb_val, None)
            pubmed_val = pubmed_dict.get(pdb_val, "")
            ph_val = ph_dict.get(pdb_val, "")
            ligand_val = ligand_dict.get(pdb_val, "")

            # Format score if present
            try:
                score_val = round(float(score_val), 3) if pd.notna(score_val) else None
            except Exception:
                score_val = None

            # Extract PEG
            cleaned_compounds, peg_val, peg_cons_val = extract_peg(comp_val)

            data.append({
                "PDB_ID": pdb_val,
                "SCORE": score_val,
                "PUBMED_ID": pubmed_val,
                "PH": ph_val,
                "LIGAND": ligand_val,
                "COMPOUNDS": cleaned_compounds,
                "PEGS": peg_val,
                "PEG_CONS": peg_cons_val
            })
            found_pdb_ids.add(pdb_val)

    # -------------------------------
    # Warn about missing PDB_IDs
    # -------------------------------
    missing_pdb_ids = set(filtered_df["PDB_ID"]) - found_pdb_ids
    if missing_pdb_ids:
        print(f"Warning: {len(missing_pdb_ids)} PDB_IDs from CSV not found in pickle:")
        for pdb in sorted(missing_pdb_ids):
            print(f"  - {pdb}")

    # -------------------------------
    # Save output CSV sorted by SCORE descending
    # -------------------------------
    filtered_structures_df = pd.DataFrame(data)
    filtered_structures_df["SCORE"] = pd.to_numeric(filtered_structures_df["SCORE"], errors="coerce")
    filtered_structures_df = filtered_structures_df.sort_values(by="SCORE", ascending=False)
    filtered_structures_df.to_csv(output_csv_file, index=False)
    print(f"Filtered {len(filtered_structures_df)} structures saved to: {output_csv_file}")