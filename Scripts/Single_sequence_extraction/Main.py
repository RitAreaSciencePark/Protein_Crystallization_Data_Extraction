from PDB_searchAPI import run_pdb_search
from Plot_crys_details import run_plot

def main():
    csv_file = run_pdb_search()
    run_plot(csv_file)
    print("Pipeline completed successfully")

if __name__ == "__main__":
    main()