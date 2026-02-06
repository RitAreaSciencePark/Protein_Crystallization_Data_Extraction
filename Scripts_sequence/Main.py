import os
from PDB_searchAPI import run_pdb_search
from Plot_crys_details import run_plot
from chem_class_reagents import process_csv

def main():

    # ---------------------------
    # Create output folder
    # ---------------------------
    base_dir = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
    output_dir = os.path.join(base_dir, "output")
    os.makedirs(output_dir, exist_ok=True)

    # Run PDB search and save CSV
    pdb_csv = run_pdb_search()
    pdb_csv_out = os.path.join(output_dir, os.path.basename(pdb_csv))
    os.replace(pdb_csv, pdb_csv_out)  # move generated CSV to output folder


    # ---------------------------
    # Step 2: Parse proteins and reagents
    # ---------------------------
    parsed_csv = os.path.join(output_dir, "protein_reagents_parsed.csv")
    process_csv(pdb_csv_out, parsed_csv)

    # Plot the data
    run_plot(pdb_csv_out)
    print("Pipeline completed successfully")

if __name__ == "__main__":
    main()
