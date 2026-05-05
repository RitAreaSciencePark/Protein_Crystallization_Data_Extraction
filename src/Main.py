import os
import tempfile
import pandas as pd

from PDB_searchAPI import search_pdb_by_sequence, filter_experimental_conditions
from plot import run_plot
from extract_structures import append_compound_to_filtered_csv
from rcsb_sequence_identity import run_and_save


def main():

    # =========================================================
    # INPUT
    # =========================================================
    sequence = input("Enter sequence: ").strip()
    seq_type_name = input("Enter a descriptive sequence type name: ").strip()

    # =========================================================
    # OUTPUT FOLDER
    # =========================================================
    base_dir = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
    output_dir = os.path.join(base_dir, "output", seq_type_name)
    os.makedirs(output_dir, exist_ok=True)

    # =========================================================
    # SAVE FASTA
    # =========================================================
    fasta_path = os.path.join(output_dir, f"{seq_type_name}_sequence.fasta")
    with open(fasta_path, "w") as f:
        f.write(f">{seq_type_name}\n")
        for i in range(0, len(sequence), 60):
            f.write(sequence[i:i+60] + "\n")

    print(f"✔ FASTA saved: {fasta_path}")

    # =========================================================
    # STEP 1 — RCSB SEQUENCE SEARCH
    # =========================================================
    rcsb_csv = os.path.join(output_dir, f"{seq_type_name}_rcsb_hits.csv")

    print("\n▶ Running RCSB search...")
    rcsb_df = run_and_save(sequence, output_csv_1=rcsb_csv)

    if rcsb_df is None or rcsb_df.empty:
        print("⚠ No RCSB hits found")
    else:
        print(f"✔ RCSB CSV saved: {rcsb_csv}")

    # =========================================================
    # STEP 2 — CRYSTALLIZATION PIPELINE
    # =========================================================
    with tempfile.NamedTemporaryFile(mode="w+", delete=False, suffix=".csv") as tmp_full_csv:
        tmp_full_csv_path = tmp_full_csv.name

    search_pdb_by_sequence(
        sequence,
        output_csv=tmp_full_csv_path,
        max_workers=6
    )

    with tempfile.NamedTemporaryFile(mode="w+", delete=False, suffix=".csv") as filtered_csv:
        filtered_csv_path = filtered_csv.name

    filter_experimental_conditions(tmp_full_csv_path, filtered_csv_path)
    os.remove(tmp_full_csv_path)

    # =========================================================
    # ADD COMPOUNDS
    # =========================================================
    structures_file = os.path.join(base_dir, "Structures", "structures.pkl")

    with tempfile.NamedTemporaryFile(mode="w+", delete=False, suffix=".csv") as tmp_cryst_csv:
         cryst_csv_path = tmp_cryst_csv.name

    append_compound_to_filtered_csv(structures_file, filtered_csv_path, cryst_csv_path)

    os.remove(filtered_csv_path)
   
    # =========================================================
    # STEP 3 — MERGE RCSB + CRYSTALLIZATION DATA
    # =========================================================
    if rcsb_df is not None and not rcsb_df.empty:

        cryst_df = pd.read_csv(cryst_csv_path)


        # cleanup temp file early
        os.remove(cryst_csv_path)

        # standardize keys
        cryst_df["PDB_ID"] = cryst_df["PDB_ID"].astype(str).str.upper().str.strip()
        rcsb_df["PDB_ID"] = rcsb_df["PDB_ID"].astype(str).str.upper().str.strip()
        cryst_df = cryst_df.drop(columns=["Score"], errors="ignore")

        merged_df = pd.merge(cryst_df, rcsb_df, on="PDB_ID", how="left")
        merged_csv = os.path.join(output_dir, f"{seq_type_name}_merged_results.csv")

        column_order = ["PDB_ID", "Entity", "Score", "Seq_id", "E-value", "Resolution", "Pubmed_id", "Method", "pH", "Temp", "Ligands", "Polymer","Assembly","pdbx_pH_range",
         "pdbx_details", "Compounds(con_unit=mM)", "PEG_Id", "PEG_con"]   

        merged_df = merged_df[column_order]

        merged_df.to_csv(merged_csv, index=False)

        print(f"\n✔ Merged CSV saved: {merged_csv}")

    else:
        print("\n⚠ No RCSB data to merge")

    # =========================================================
    # STEP 4 — PLOT CRYSTALLIZATION DATA
    # =========================================================
    run_plot(merged_csv, seq_type_name)

    # =========================================================
    # DONE
    # =========================================================
    print("\n✔ Pipeline completed successfully")
    print(f"📁 Output folder: {output_dir}")


# =========================================================
# RUN
# =========================================================
if __name__ == "__main__":
    main()