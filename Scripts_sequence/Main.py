import os
import pandas as pd
from PDB_searchAPI import search_pdb_by_sequence, filter_experimental_conditions
from plot import run_plot
from extract_structures import extract_structures_with_metadata

# -------------------------------
# Main pipeline
# -------------------------------
def main():
    sequence = input("Enter sequence: ").strip()
    seq_type = input("Enter sequence type (protein/dna/rna): ").strip().lower()
    seq_type_name = input("Enter a descriptive sequence type name: ").strip()

    # Create output folder
    base_dir = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
    output_dir = os.path.join(base_dir, "output", seq_type_name)
    os.makedirs(output_dir, exist_ok=True)

    # Save FASTA
    fasta_path = os.path.join(output_dir, f"{seq_type_name}_sequence.fasta")
    with open(fasta_path, "w") as fasta_file:
        fasta_file.write(f">{seq_type_name}\n")
        for i in range(0, len(sequence), 60):
            fasta_file.write(sequence[i:i+60] + "\n")
    print(f"✔ Input sequence saved as FASTA: {fasta_path}")

    # Run PDB search
    full_csv = os.path.join(output_dir, f"{seq_type_name}_pdb_mmcif_full.csv")
    full_csv_result = search_pdb_by_sequence(sequence, seq_type, output_csv=full_csv, keep_all=True)

    if full_csv_result is None:
        print(f"⚠ No PDB entries found for this {seq_type} sequence. Pipeline stopped.")
        return

    # Filter experimental conditions
    filtered_csv = os.path.join(output_dir, f"{seq_type_name}_pdb_mmcif_filtered.csv")
    filter_experimental_conditions(full_csv, filtered_csv)

    # -------------------------------
    # Extract structures with metadata from pickle
    # -------------------------------
    structures_pickle_file = os.path.join(base_dir, "Structures", "structures.pkl")
    output_structures_csv = os.path.join(output_dir, f"{seq_type_name}_reagents.csv")

    extract_structures_with_metadata(structures_pickle_file, filtered_csv, output_structures_csv)

    # -------------------------------
    # Parse reagents (if needed)
    # -------------------------------
    # parsed_csv = os.path.join(output_dir, f"{seq_type_name}_reagents_parsed.csv")
    # json_path = os.path.join(base_dir, "Data", "reagents.json")  # Ensure your updated JSON is here
    # process_csv(filtered_csv, parsed_csv, json_path)

    # Generate plots
    run_plot(filtered_csv)

    print(f"\n✔ Pipeline completed successfully.")
    print(f"✔ All outputs saved in {output_dir}")

# -------------------------------
# Run main
# -------------------------------
if __name__ == "__main__":
    main()