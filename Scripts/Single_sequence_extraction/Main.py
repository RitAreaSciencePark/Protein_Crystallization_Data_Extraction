from PDB_searchAPI import run_pdb_search
from Plot_crys_details import run_plot
from chem_class_reagents import process_csv

def main():
    # Run PDB search and save CSV
    pdb_csv = run_pdb_search()

    # Parse proteins and reagents into final CSV
    processed_csv = process_csv(pdb_csv, "protein_reagents_parsed.csv")

    # Plot the data
    run_plot(pdb_csv)
    print("Pipeline completed successfully")

if __name__ == "__main__":
    main()
