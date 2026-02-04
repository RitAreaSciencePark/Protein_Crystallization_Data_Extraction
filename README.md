
# Protein Crystallization Data Extraction (PCDE)

The protein crystallization data extraction (PCDE) is a python based tool that automatically retrieves, organizes, and analyses data related to protein crystallization conditions uploaded in the Protein Data Bank (PDB). The tool extract all the data related to experimental crystallization conditions of all identified homologous proteins to a single sequence or multiple sequences stored in a fasta file used as input to the codes. 

The scripts the script folder are divided into two main folders: the single sequence extraction (Single_sequence_extraction folder) and the multiple sequence extraction (Scripts_fasta_file folder)

## Scripts_fasta_file folder

The scripts_fasta_file folder contains the scripts used to extract the data from the PDB and save it to a csv file. The scripts are:

### main.py

This script is the main script that calls import *read_and_clean_fasta_file*, *write_fasta_dict* functions from *cleaning_and_read_fasta_file.py* and *fasta_to_pdb_csv* from *extract_data_fasta.py*.

- cleaning_and_read_fasta_file.py
   This script aims to: read the raw fasta file, remove all '-' from sequences, and returns a dictionary mapping sequence id to sequence. It handles multiple sequences stored in a single fasta file.

- extract_data_fasta.py
   The cleaned fasta file is then used at this stage for a RCSB PDB search to find the sequence identity matching and extract the crystallization data from the PDB. The extracted data is then saved to a csv file.
   
   The script uses the RSCB search API to search for the sequence identity matching through a JSON query based on the sequence and filter the results to only include the structures determined by X-ray diffraction. In order to have the maximum number of structures the identity_cutoff is set to 30%.

   For each PDB hit, it get the PDB ID and search the score.

   The function '''extract_crystallization(pdb_id)'''pulls all experimental data from the mmCIF file by parsing the file using *gemmi*, accessing the main block with *doc.sole_block()*, and then returns a dictionary with the experimental details

   The csv file gives the query_ID, pdb_id, score, and the experimental crystallization data of similar sequences.

to run the script: '''python main.py <input_fasta.fasta>'''

### main2.py

This script is the main script that calls import *read_and_clean_fasta_file*, *write_fasta_dict* functions from *cleaning_and_read_fasta_file.py* and *fasta_to_pdb_csv* from *fasta_to_pdb*.

The main difference with main.py is that the csv file contains the combined data of all the sequences stored in the fasta file without repetition of PDB ID.

to run the script: '''python main2.py <input_fasta.fasta>'''


## single_sequence_extraction folder

The main.py script handles the PDB-search with a single sequence as input, the extraction of protein/reagent and their concentration from the pdbx_details and plot the pH against the temperature, the method (sitting/hanging drop).

The script imports *run_pdb_search*, *run_plot* and *process_csv* functions from *PDB_searchAPI.py*, *plot_crys_details.py* and *chem_class_reagents.py*.

   - The PDB_searchAPI.py is the same as in the multiple sequence extraction. the results are saved to a csv file called "pdb_mmcif_extracted.csv".

   - The plot_crys_details.py
      This script load the "pdb_mmcif_extracted.csv" file from the pdb search,
      filter for crystallization methods ("hanging", or "sitting" drops),
      Parse the pH - range: if the pH is 6.0 - 8.0, it parses it into two columns: pH_low and pH_high and compute the errors,
      Uses *viridis* colormap to represent the score and normalize the scores to map them,
      output the plot as pdf file.
      

   - The chem_class_reagents.py 
       Create a reagent and protein dictionarie,
       Handles multiples proteins in one pdbx_details entry,
       Extracts reagents and their concentrations from the pdbx_details entry,
       Classify reagents into buffer, salt, peg, ligand and solvent,
       Normalizes chemical synonyms and formula,
       output a csv file with the reagents and their concentrations.

To run the script: '''python main.py'''  then enter the sequence of a protein of interest.

It will create an output folder where the csv file with the reagents and their concentrations, the pdf plot and the pdb_mmcif_extracted.csv file will be saved.
       
