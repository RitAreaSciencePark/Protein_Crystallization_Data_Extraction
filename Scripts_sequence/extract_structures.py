import pickle
import pandas as pd
import types
import sys

def append_compound_to_filtered_csv(structures_file, filtered_csv_file, output_csv_file):
    """
    Append COMPOUND column from pickle to filtered CSV based on matching PDB_ID.
    Saves output to output_csv_file (can overwrite filtered_csv_file).
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

    if not structures_list:
        print("⚠ Pickle file is empty. No COMPOUND info to append.")
        return

    # -------------------------------
    # Load filtered CSV
    # -------------------------------
    filtered_df = pd.read_csv(filtered_csv_file, on_bad_lines="skip")
    filtered_df["PDB_ID"] = filtered_df["PDB_ID"].astype(str).str.strip().str.upper()

    # -------------------------------
    # Detect PDB_ID and COMPOUND attributes in pickle
    # -------------------------------
    first_obj = structures_list[0]
    keys = first_obj.__dict__.keys()

    pdb_attr = next((k for k in keys if "pdb" in k.lower()), None)
    comp_attr = next((k for k in keys if "compound" in k.lower()), None)

    if not pdb_attr or not comp_attr:
        print("⚠ Pickle does not contain PDB_ID or COMPOUND attributes. Skipping append.")
        return

    # Build lookup dictionary: PDB_ID -> COMPOUND
    compound_dict = {
        str(getattr(obj, pdb_attr, "")).strip().upper(): getattr(obj, comp_attr, "")
        for obj in structures_list
    }

    # Map COMPOUND values to filtered CSV
    filtered_df["CRYSTALLIZATION_COCKTAILS"] = filtered_df["PDB_ID"].map(compound_dict).fillna("")

    # -------------------------------
    # Save CSV
    # -------------------------------
    filtered_df.to_csv(output_csv_file, index=False)
    print(f"✔ COMPOUND column appended. CSV saved to: {output_csv_file}")