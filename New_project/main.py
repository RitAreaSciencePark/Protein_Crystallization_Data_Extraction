#!/usr/bin/env python3

"""
main.py
=======

Orchestrates the two-stage pipeline:

  1. pdb_sequence_search.py
     Given a protein/nucleic-acid sequence, search RCSB PDB for X-ray
     homologs, pull crystallization metadata for each hit, and filter down
     to entries with real experimental conditions -> writes Output.csv
     (the only file that stage produces).

  2. compound_extraction.py
     Reads Output.csv, extracts (compound, amount, unit) triples from its
     `pdbx_details` column (dictionary/regex pass, optionally with an LLM
     fallback that also *learns* new compounds for next time), and writes
     the expanded compound-level CSV.

  3. plot.py
     Reads Output_compounds.csv and generates the pH-vs-Temperature plot,
     pH-vs-PEG-concentration plot, and a colored PDF summary table, plus
     Grouped_conditions.csv (one row per unique condition, with a row_id
     shared between plot points and table rows -- this is the file the
     Django/Plotly web app reads to implement click-a-point-to-jump-to-
     its-table-row, since static PNG/PDF can't do that on their own).

Every run asks for (or takes via --protein-name) a protein name, and all
output files for that run are saved inside a folder named after it, e.g.:

    Lysozyme/
      Output.csv
      Output_compounds.csv

Usage
-----
    # interactive prompt for both the protein name and the sequence
    python main.py

    # everything given directly
    python main.py --protein-name Lysozyme --sequence "MKTAYIAKQ..."

    # from a FASTA file, with the LLM fallback enabled for compound extraction
    python main.py --protein-name Lysozyme --fasta my_protein.fasta --llm-fallback

Both stages' options are exposed here; run `python main.py --help` for the
full list.
"""

import argparse
import os
import re
import sys
from pdb_sequence_search import (find_homologs_with_conditions, filter_experimental_conditions, prompt_for_sequence, read_fasta)
from compound_extraction import process_csv
from plot import run_plot

def prompt_for_protein_name() -> str:
    """Interactively ask for the protein name from the console."""
    name = input("Enter the protein name (used to name the output folder): ").strip()
    return name

def sanitize_folder_name(name: str) -> str:
    """Turn a free-text protein name into a safe, filesystem-friendly
    folder name: spaces -> underscores, anything else non-alphanumeric
    (besides '_' and '-') stripped out."""
    name = name.strip().replace(" ", "_")
    name = re.sub(r"[^A-Za-z0-9_\-]", "", name)
    return name or "protein_output"  # fallback if sanitizing empties the string

def main():
    p = argparse.ArgumentParser( description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)

    # --- protein name / output folder ---
    p.add_argument("--protein-name", default=None,
                   help="Name of the protein being searched. Used to name the output "
                        "folder. If omitted, you'll be prompted for it interactively.")
    
    # --- Save Fasta
    fasta_path = os.path.join(output_dir, f"{protein_name}_sequence.fasta")
    with open(fasta_path, "w") as f:
        f.write(f">{protein_name}\n")
        for i in range(0, len(sequence), 60):
            f.write(sequence[i:i+60] + "\n")

    print(f"✔ FASTA saved: {fasta_path}")

    # --- sequence input (stage 1) ---
    seq_group = p.add_mutually_exclusive_group(required=False)
    seq_group.add_argument("--sequence", help="Protein/nucleic acid sequence, as a raw string.")
    seq_group.add_argument("--fasta", help="Path to a FASTA file containing the query sequence.")
    p.add_argument("--sequence-type", default="protein", choices=["protein", "dna", "rna"],
                   help="Sequence type for the search (default: protein).")
    p.add_argument("--identity", type=float, default=0.5,
                   help="Minimum fractional sequence identity, 0-1 (default: 0.5).")
    p.add_argument("--evalue", type=float, default=1.0,
                   help="Maximum E-value cutoff (default: 1.0).")
    p.add_argument("--max-hits", type=int, default=100,
                   help="Maximum number of PDB entries to retrieve (default: 100).")
    p.add_argument("--sleep", type=float, default=0.2,
                   help="Seconds to sleep between per-entry API calls (default: 0.2).")
    p.add_argument("--quiet", action="store_true", help="Suppress progress messages.")

    # --- stage 1 output / stage 2 input (filenames only -- always placed
    # inside the per-protein output folder, see below) ---
    p.add_argument("--out", default="Output.csv",
                   help="Filename for the filtered crystallization-conditions CSV "
                        "(default: Output.csv), saved inside the protein's output "
                        "folder. This is stage 1's output AND stage 2's input.")

    # --- compound extraction (stage 2) ---
    p.add_argument("--compounds-out", default="Output_compounds.csv",
                   help="Filename for the final compound-level CSV "
                        "(default: Output_compounds.csv), saved inside the protein's "
                        "output folder.")
    p.add_argument("--llm-fallback", action="store_true",
                   help="Use Claude to resolve pdbx_details clauses the dictionary/regex "
                        "pass can't match, and learn any new compounds it finds.")
    p.add_argument("--api-key", default=None,
                   help="Anthropic API key for --llm-fallback (or set ANTHROPIC_API_KEY).")

    # --- allow skipping stage 1 if you already have an Output.csv ---
    p.add_argument("--skip-search", action="store_true",
                   help="Skip stage 1 (sequence search) and run compound extraction "
                        "directly on an existing Output.csv inside the protein's "
                        "output folder.")
    p.add_argument("--skip-plots", action="store_true",
                   help="Skip stage 3 (plots/PDF table generation).")

    args = p.parse_args()

    # ---------------- Protein name -> output folder ----------------
    protein_name = args.protein_name or prompt_for_protein_name()
    protein_name = protein_name.strip()
    if not protein_name:
        p.error("No protein name provided.")

    output_dir = sanitize_folder_name(protein_name)
    os.makedirs(output_dir, exist_ok=True)
    print(f"=== Output folder: {os.path.abspath(output_dir)} ===", file=sys.stderr)

    out_path = os.path.join(output_dir, args.out)
    compounds_out_path = os.path.join(output_dir, args.compounds_out)

    if not args.skip_search:
        # ---------------- Stage 1: sequence search -> Output.csv ----------------
        if args.sequence:
            sequence = args.sequence
        elif args.fasta:
            sequence = read_fasta(args.fasta)
        else:
            sequence = prompt_for_sequence()

        sequence = sequence.strip().upper()
        if not sequence:
            p.error("No sequence provided (nothing entered / FASTA file was empty).")

        print("=== Stage 1: searching RCSB PDB for homologs ===", file=sys.stderr)
        rows = find_homologs_with_conditions(sequence, identity_cutoff=args.identity, evalue_cutoff=args.evalue, sequence_type=args.sequence_type,
            max_hits=args.max_hits, sleep=args.sleep, verbose=not args.quiet,)
        # rows stays in memory -- filter_experimental_conditions() is the only
        # function that touches disk here, and out_path is the only file it writes.
        filter_experimental_conditions(rows, output_csv=out_path)
    else:
        print(f"=== Stage 1 skipped: using existing {out_path} ===", file=sys.stderr)

    # ---------------- Stage 2: Output.csv -> compound-level CSV ----------------
    print("=== Stage 2: extracting compounds/concentrations ===", file=sys.stderr)
    process_csv(
        out_path,                    # <-- Output.csv from stage 1 is stage 2's input
        compounds_out_path,
        use_llm_fallback=args.llm_fallback,
        api_key=args.api_key,
    )
    print(f"Done. Crystallization conditions: {out_path}", file=sys.stderr)
    print(f"Done. Compound-level results: {compounds_out_path}", file=sys.stderr)

    # ---------------- Stage 3: Output_compounds.csv -> plots + PDF table ----------------
    if not args.skip_plots:
        print("=== Stage 3: generating plots and PDF table ===", file=sys.stderr)
        plot_paths = run_plot(output_dir, protein_name, compounds_csv_name=args.compounds_out)
        for label, path in plot_paths.items():
            print(f"Done. {label}: {path}", file=sys.stderr)
    else:
        print("=== Stage 3 skipped ===", file=sys.stderr)


if __name__ == "__main__":
    main()