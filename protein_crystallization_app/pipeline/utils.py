import os
import tempfile
from .PDB_searchAPI import search_pdb_by_sequence, filter_experimental_conditions
from .plot import run_plot
from .extract_structures import append_compound_to_filtered_csv

def run_pipeline(sequence, seq_type_name, media_root):
    """
    Runs the protein crystallization pipeline and saves outputs directly into media_root/<seq_type_name>/
    Returns a dictionary of file paths for CSV, FASTA, PNGs, and PDFs
    """
    # -------------------------------
    # Output folder under MEDIA_ROOT
    # -------------------------------
    output_dir = os.path.join(media_root, seq_type_name)
    os.makedirs(output_dir, exist_ok=True)

    # -------------------------------
    # Save input sequence as FASTA
    # -------------------------------
    fasta_path = os.path.join(output_dir, f"{seq_type_name}_sequence.fasta")
    with open(fasta_path, "w") as fasta_file:
        fasta_file.write(f">{seq_type_name}\n")
        for i in range(0, len(sequence), 60):
            fasta_file.write(sequence[i:i+60] + "\n")

    # -------------------------------
    # Run PDB search and temporary CSVs
    # -------------------------------
    with tempfile.NamedTemporaryFile(mode="w+", delete=False, suffix=".csv") as tmp_full_csv:
        tmp_full_csv_path = tmp_full_csv.name

    full_csv_result = search_pdb_by_sequence(sequence, output_csv=tmp_full_csv_path, keep_all=True)
    if full_csv_result is None:
        os.remove(tmp_full_csv_path)
        return None  # No PDB entries

    with tempfile.NamedTemporaryFile(mode="w+", delete=False, suffix=".csv") as filtered_csv:
        filtered_csv_path = filtered_csv.name

    filter_experimental_conditions(tmp_full_csv_path, filtered_csv_path)
    os.remove(tmp_full_csv_path)

    # -------------------------------
    # Append COMPOUND column and save final CSV
    # -------------------------------
    structures_file = os.path.join(os.path.dirname(os.path.dirname(os.path.abspath(__file__))), "Structures", "structures.pkl")
    output_csv_file = os.path.join(output_dir, f"{seq_type_name}_crystallization_data.csv")
    append_compound_to_filtered_csv(structures_file, filtered_csv_path, output_csv_file)
    os.remove(filtered_csv_path)

    # -------------------------------
    # Generate plots
    # -------------------------------
    run_plot(output_csv_file)
    # Assuming run_plot saves a PNG and PDF with consistent names:
    full_plot_file = os.path.join(output_dir, f"{seq_type_name}_FULL_PLOT.png")
    first10_pdf_file = os.path.join(output_dir, f"{seq_type_name}_FIRST10_TABLE.pdf")

    return {
        "csv": output_csv_file,
        "fasta": fasta_path,
        "pngs": [full_plot_file],
        "pdfs": [first10_pdf_file]
    }