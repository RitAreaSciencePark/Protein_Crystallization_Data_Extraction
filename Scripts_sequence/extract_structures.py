# extract_structures.py
import sys
import types
import pickle
import pandas as pd
from pathlib import Path

def extract_structures_with_metadata(structures_file, filtered_csv_file, output_csv_file):
    """
    Extract PDB_ID, COMPOUNDS from pickle and merge SCORE, PH, PUBMED_ID from filtered CSV.
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

    if "PDB_ID" not in filtered_df.columns:
        raise ValueError("Filtered CSV must have a 'PDB_ID' column.")

    # Correct: use .str.upper()
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

    # Build lookup dictionaries, handle missing columns
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
    if not pdb_attr:
        raise ValueError("No attribute resembling PDB_ID found in pickle objects.")
    comp_attr = next((k for k in keys if 'compound' in k.lower()), None)
    if not comp_attr:
        raise ValueError("No attribute resembling COMPOUNDS found in pickle objects.")

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

            data.append({
                "PDB_ID": pdb_val,
                "SCORE": score_val,
                "PUBMED_ID": pubmed_val,
                "PH": ph_val,
                "LIGAND": ligand_val,
                "COMPOUNDS": comp_val
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