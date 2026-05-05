def run_pipeline(sequence, seq_type_name, base_output_dir):
    import os
    import tempfile
    import pandas as pd
    from pdf2image import convert_from_path
    from .PDB_searchAPI import search_pdb_by_sequence, filter_experimental_conditions
    from .plot import run_plot
    from .extract_structures import append_compound_to_filtered_csv
    from .rcsb_sequence_identity import run_and_save
    import shutil

    output_dir = os.path.join(base_output_dir, seq_type_name)

    # 🔥 delete previous results
    if os.path.exists(output_dir):
        shutil.rmtree(output_dir)

    # recreate fresh folder
    os.makedirs(output_dir, exist_ok=True)

    output_dir = os.path.join(base_output_dir, seq_type_name)
    os.makedirs(output_dir, exist_ok=True)

    # FASTA
    fasta_path = os.path.join(output_dir, f"{seq_type_name}_sequence.fasta")
    with open(fasta_path, "w") as f:
        f.write(f">{seq_type_name}\n")
        for i in range(0, len(sequence), 60):
            f.write(sequence[i:i+60] + "\n")

    # STEP 1 — RCSB
    rcsb_csv = os.path.join(output_dir, f"{seq_type_name}_rcsb_hits.csv")
    rcsb_df = run_and_save(sequence, output_csv_1=rcsb_csv)

    if rcsb_df is None or rcsb_df.empty:
        return None

    # STEP 2 — crystallization
    with tempfile.NamedTemporaryFile(delete=False, suffix=".csv") as tmp_full:
        tmp_full_path = tmp_full.name

    search_pdb_by_sequence(sequence, output_csv=tmp_full_path, max_workers=6)

    with tempfile.NamedTemporaryFile(delete=False, suffix=".csv") as tmp_filtered:
        filtered_path = tmp_filtered.name

    filter_experimental_conditions(tmp_full_path, filtered_path)
    os.remove(tmp_full_path)

    structures_file = os.path.join(os.path.dirname(__file__), "Structures", "structures.pkl")

    with tempfile.NamedTemporaryFile(delete=False, suffix=".csv") as tmp_cryst:
        cryst_path = tmp_cryst.name

    append_compound_to_filtered_csv(structures_file, filtered_path, cryst_path)
    os.remove(filtered_path)

    # STEP 3 — merge
    cryst_df = pd.read_csv(cryst_path)
    os.remove(cryst_path)

    cryst_df["PDB_ID"] = cryst_df["PDB_ID"].astype(str).str.upper().str.strip()
    rcsb_df["PDB_ID"] = rcsb_df["PDB_ID"].astype(str).str.upper().str.strip()
    cryst_df = cryst_df.drop(columns=["Score"], errors="ignore")

    merged_df = pd.merge(cryst_df, rcsb_df, on="PDB_ID", how="left")

    merged_csv = os.path.join(output_dir, f"{seq_type_name}_merged_results.csv")

    column_order = [
        "PDB_ID", "Entity", "Score", "Seq_id", "E-value",
        "Resolution", "Pubmed_id", "Method", "pH", "Temp",
        "Ligands", "Polymer", "Assembly", "pdbx_pH_range",
        "pdbx_details", "Compounds(con_unit=mM)", "PEG_Id", "PEG_con"
    ]

    merged_df = merged_df[column_order]
    merged_df.to_csv(merged_csv, index=False)

    # STEP 4 — plots
    run_plot(merged_csv, seq_type_name)

    # collect generated files
    png_files = [os.path.join(output_dir, f) for f in os.listdir(output_dir) if f.endswith(".png")]
    pdf_files = [os.path.join(output_dir, f) for f in os.listdir(output_dir) if f.endswith(".pdf")]
    
    return {
    "merged_csv": merged_csv,
    "fasta": fasta_path,
    "pngs": png_files,   
    "pdfs" : pdf_files          
}
