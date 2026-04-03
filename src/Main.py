import os
from PDB_searchAPI import search_pdb_by_sequence, filter_experimental_conditions
from plot import run_plot
from extract_structures import append_compound_to_filtered_csv  # COMPOUND + Excel function
import tempfile

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

    full_csv_result = search_pdb_by_sequence(sequence, output_csv=tmp_full_csv_path, max_workers=6)

    if full_csv_result is None:
        print(f"⚠ No PDB entries found for this {seq_type_name} sequence. Pipeline stopped.")
        os.remove(tmp_full_csv_path)
        return

    # -------------------------------
    # Filter experimental conditions (into temporary file)
    # -------------------------------
    with tempfile.NamedTemporaryFile(mode="w+", delete=False, suffix=".csv") as filtered_csv:
        filtered_csv_path = filtered_csv.name

    filter_experimental_conditions(tmp_full_csv_path, filtered_csv_path)
    os.remove(tmp_full_csv_path)  # remove temporary full CSV

    # -------------------------------
    # Append COMPOUND column and save only final CSV & Excel
    # -------------------------------
    structures_file = os.path.join(base_dir, "Structures", "structures.pkl")

   
    # Save COMPOUND-augmented CSV 
    output_csv_file = os.path.join(output_dir, f"{seq_type_name}_crystallization_data.csv")
 
    append_compound_to_filtered_csv(structures_file, filtered_csv_path, output_csv_file)
    os.remove(filtered_csv_path)  # Clean up temporary filtered CSV

    # -------------------------------
    # Generate plots using final CSV
    # ------------------------------- 
    run_plot(output_csv_file)

    print(f"\n✔ Pipeline completed successfully.")
   
# -------------------------------
# Run main
# -------------------------------
if __name__ == "__main__":
    main()