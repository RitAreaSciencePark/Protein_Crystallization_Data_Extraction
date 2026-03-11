import os
from PDB_searchAPI import search_pdb_by_sequence, filter_experimental_conditions
from plot import run_plot
from extract_structures import append_compound_to_filtered_csv  # COMPOUND-only function
import tempfile

# -------------------------------
# Main pipeline
# -------------------------------
def main():
    # -------------------------------
    # Input
    # -------------------------------
    sequence = input("Enter sequence: ").strip()
    seq_type_name = input("Enter a descriptive sequence type name: ").strip()

    # -------------------------------
    # Create output folder
    # -------------------------------
    base_dir = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
    output_dir = os.path.join(base_dir, "output", seq_type_name)
    os.makedirs(output_dir, exist_ok=True)

    # -------------------------------
    # Save input sequence as FASTA
    # -------------------------------
    fasta_path = os.path.join(output_dir, f"{seq_type_name}_sequence.fasta")
    with open(fasta_path, "w") as fasta_file:
        fasta_file.write(f">{seq_type_name}\n")
        for i in range(0, len(sequence), 60):
            fasta_file.write(sequence[i:i+60] + "\n")
    print(f"✔ Input sequence saved as FASTA: {fasta_path}")

    # -------------------------------
    # Run PDB search and write to a temporary full CSV
    # -------------------------------
    with tempfile.NamedTemporaryFile(mode="w+", delete=False, suffix=".csv") as tmp_full_csv:
        tmp_full_csv_path = tmp_full_csv.name

    full_csv_result = search_pdb_by_sequence(sequence, output_csv=tmp_full_csv_path, keep_all=True)

    if full_csv_result is None:
        print(f"⚠ No PDB entries found for this {seq_type_name} sequence. Pipeline stopped.")
        os.remove(tmp_full_csv_path)
        return

    # -------------------------------
    # Filter experimental conditions
    # -------------------------------
    filtered_csv = os.path.join(output_dir, f"{seq_type_name}_pdb_mmcif_filtered.csv")
    filter_experimental_conditions(tmp_full_csv_path, filtered_csv)
    os.remove(tmp_full_csv_path)  # Clean up temporary CSV

    # -------------------------------
    # Append COMPOUND column from pickle
    # -------------------------------
    structures_pickle_file = os.path.join(base_dir, "Structures", "structures.pkl")
    append_compound_to_filtered_csv(structures_pickle_file, filtered_csv, filtered_csv)  # overwrite same file

    # -------------------------------
    # Generate plots using final CSV
    # -------------------------------
    run_plot(filtered_csv)

    print(f"\n✔ Pipeline completed successfully.")
    print(f"✔ All outputs saved in {output_dir}")


# -------------------------------
# Run main
# -------------------------------
if __name__ == "__main__":
    main()