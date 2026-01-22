# main.py

import sys
import csv

from cleaning_and_read_fasta_file import read_and_clean_fasta_file, write_fasta_dict
from Reading_fasta import fasta_to_pdb_csv


def main():
    # Check if a FASTA file was provided
    if len(sys.argv) < 2:
        print("Usage: python main.py <input_fasta>")
        sys.exit(1)

    input_fasta = sys.argv[1]  # <-- This reads the file path from terminal
   
    # Read and clean sequences
    fasta_dict = read_and_clean_fasta_file(input_fasta)


    # Print cleaned sequences to terminal
    #for query_id, sequence in fasta_dict.items():
       # print(f">{query_id}")
       # print(sequence)

    # Optional: save to a new file
    cleaned_fasta = "cleaned.fasta"
    write_fasta_dict(fasta_dict, cleaned_fasta)
    print(f"\nCleaned sequences saved to {cleaned_fasta}")

    # Run PDB search using cleaned sequences (NO re-reading FASTA)
    output_csv = "pdb2_results.csv" # single CSV output
    fasta_to_pdb_csv(
        fasta_file=cleaned_fasta,
        output_csv=output_csv,
        identity_cutoff=0.2,
        rows=1000
    )

if __name__ == "__main__":
    main()
