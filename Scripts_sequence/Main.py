import os
from PDB_searchAPI import run_pdb_search, filter_experimental_conditions
from plot import run_plot
from chem_class_reagents import process_csv

def main():
    protein_name = input("Enter a name for this protein (for folder and filenames):\n").strip()
    if not protein_name:
        protein_name = "protein"
        
    sequence = input("Enter protein sequence (single-letter amino acid code, no FASTA header):\n").strip().upper()
    if not sequence:
        raise ValueError("Sequence cannot be empty")

    # Create protein output folder
    base_dir = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
    protein_output_dir = os.path.join(base_dir, "output", protein_name)
    os.makedirs(protein_output_dir, exist_ok=True)

    # Save input sequence as FASTA
    fasta_path = os.path.join(protein_output_dir, f"{protein_name}_sequence.fasta")
    with open(fasta_path, "w") as fasta_file:
        fasta_file.write(f">{protein_name}\n")
        for i in range(0, len(sequence), 60):
            fasta_file.write(sequence[i:i+60] + "\n")
    print(f"✔ Input sequence saved as FASTA: {fasta_path}")

    # Run PDB search (full CSV)
    full_csv = os.path.join(protein_output_dir, f"{protein_name}_pdb_mmcif_full.csv")
    run_pdb_search(sequence, output_csv=full_csv, keep_all=True)

    # Filter CSV
    filtered_csv = os.path.join(protein_output_dir, f"{protein_name}_pdb_mmcif_filtered.csv")
    filter_experimental_conditions(full_csv, filtered_csv)

    # Parse and plot
    parsed_csv = os.path.join(protein_output_dir, f"{protein_name}_protein_reagents_parsed.csv")
    process_csv(filtered_csv, parsed_csv)
    run_plot(filtered_csv)

    print(f"✔ Pipeline completed successfully. All outputs saved in {protein_output_dir}")

if __name__ == "__main__":
    main()
