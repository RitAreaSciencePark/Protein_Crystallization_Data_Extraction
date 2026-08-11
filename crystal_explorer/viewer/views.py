import os
import re
import sys
import traceback

from django.conf import settings
from django.http import FileResponse, Http404
from django.shortcuts import render, redirect
from django.urls import reverse

from .forms import SequenceSearchForm
from .data import (
    build_peg_plot, build_temp_plot, build_table_rows, build_method_legend,
    load_conditions, score_gradient_stops, TABLE_COLUMNS, SCORE_MIN, SCORE_MAX, SCORE_TICKS,
)

# Make the pipeline package importable (pipeline/ sits at the project root).
sys.path.insert(0, str(settings.BASE_DIR))
from pipeline.pdb_sequence_search import find_homologs_with_conditions, filter_experimental_conditions
from pipeline.compound_extraction import process_csv
from pipeline.plot import run_plot


def sanitize_folder_name(name: str) -> str:
    """Same sanitizer as the CLI's main.py, kept in sync here."""
    name = name.strip().replace(" ", "_")
    name = re.sub(r"[^A-Za-z0-9_\-]", "", name)
    return name or "protein_output"


def index(request):
    """GET: show the search form. POST: run the full pipeline (save FASTA
    -> search -> filter -> compound extraction -> plots/table)
    synchronously, then redirect to the results page for that protein."""
    if request.method == "POST":
        form = SequenceSearchForm(request.POST)
        if form.is_valid():
            protein_name = form.cleaned_data["protein_name"]
            folder_name = sanitize_folder_name(protein_name)
            output_dir = os.path.join(settings.PIPELINE_OUTPUT_DIR, folder_name)
            os.makedirs(output_dir, exist_ok=True)

            try:
                sequence = form.cleaned_data["sequence"]

                # Save the input sequence as a FASTA file in the protein's
                # output folder -- this is also offered as a static export
                # on the results page.
                fasta_path = os.path.join(output_dir, f"{folder_name}.fasta")
                with open(fasta_path, "w") as f:
                    f.write(f">{folder_name}\n{sequence}\n")

                rows = find_homologs_with_conditions(
                    sequence,
                    identity_cutoff=form.cleaned_data["identity"],
                    evalue_cutoff=form.cleaned_data["evalue"],
                    sequence_type=form.cleaned_data["sequence_type"],
                    max_hits=form.cleaned_data["max_hits"],
                    verbose=True,
                )
                out_path = os.path.join(output_dir, "Output.csv")
                filter_experimental_conditions(rows, output_csv=out_path)

                compounds_out_path = os.path.join(output_dir, "Output_compounds.csv")
                process_csv(
                    out_path, compounds_out_path,
                    use_llm_fallback=form.cleaned_data["llm_fallback"],
                )

                run_plot(output_dir, folder_name)

            except Exception as e:
                return render(request, "viewer/explorer.html", {
                    "form": form,
                    "run_error": f"Pipeline failed: {e}",
                    "run_error_detail": traceback.format_exc(),
                })

            return redirect(reverse("viewer:results", args=[folder_name]))
    else:
        form = SequenceSearchForm()

    return render(request, "viewer/explorer.html", {"form": form})


def results(request, protein_name):
    """Show the plots + table for a protein that's already been run."""
    output_dir = os.path.join(settings.PIPELINE_OUTPUT_DIR, protein_name)
    grouped_csv_path = os.path.join(output_dir, "Grouped_conditions.csv")

    error = None
    peg_plot_div = None
    temp_plot_div = None
    table_rows = []
    method_legend = []

    if not os.path.exists(grouped_csv_path):
        error = (
            f"No results found for '{protein_name}' yet "
            f"(expected {grouped_csv_path}). Run a search first."
        )
    else:
        try:
            import plotly.offline as pyo
            df = load_conditions(grouped_csv_path)

            # Both plots come from the same df, so one symbol assignment
            # covers both -- see build_method_legend / _method_symbol_map.
            method_legend = build_method_legend(df)

            plot_config = {"displaylogo": False, "responsive": True}

            peg_fig = build_peg_plot(df)
            peg_plot_div = pyo.plot(
                peg_fig, output_type="div", include_plotlyjs=False, config=plot_config
            )

            temp_fig = build_temp_plot(df)
            temp_plot_div = pyo.plot(
                temp_fig, output_type="div", include_plotlyjs=False, config=plot_config
            )

            table_rows = build_table_rows(df)
        except Exception as e:
            error = f"Could not load or process {grouped_csv_path}: {e}"

    pdf_path = os.path.join(output_dir, f"{protein_name}_Cryst_cocktail_Table.pdf")
    peg_png_path = os.path.join(output_dir, f"{protein_name}_PEG.png")
    temp_png_path = os.path.join(output_dir, f"{protein_name}_TEMP.png")
    fasta_path = os.path.join(output_dir, f"{protein_name}.fasta")
    compounds_csv_path = os.path.join(output_dir, "Output_compounds.csv")

    return render(request, "viewer/explorer.html", {
        "form": SequenceSearchForm(initial={"protein_name": protein_name}),
        "protein_name": protein_name,
        "peg_plot_div": peg_plot_div,
        "temp_plot_div": temp_plot_div,
        "method_legend": method_legend,
        "score_gradient_stops": score_gradient_stops(),
        "score_min": SCORE_MIN,
        "score_max": SCORE_MAX,
        "score_ticks": SCORE_TICKS,
        "table_rows": table_rows,
        "table_columns": TABLE_COLUMNS,
        "error": error,
        "csv_path": grouped_csv_path,
        "has_pdf": os.path.exists(pdf_path),
        "has_peg_png": os.path.exists(peg_png_path),
        "has_temp_png": os.path.exists(temp_png_path),
        "has_fasta": os.path.exists(fasta_path),
        "has_compounds_csv": os.path.exists(compounds_csv_path),
    })


_DOWNLOAD_KINDS = {
    "pdf": (lambda name: f"{name}_Cryst_cocktail_Table.pdf", "application/pdf"),
    "peg": (lambda name: f"{name}_PEG.png", "image/png"),
    "temp": (lambda name: f"{name}_TEMP.png", "image/png"),
    "fasta": (lambda name: f"{name}.fasta", "text/plain"),
    "compounds": (lambda name: "Output_compounds.csv", "text/csv"),
}


def download(request, protein_name, kind):
    """Serve one of the static exports (PDF table, PEG/temp plots, input
    FASTA, compounds CSV) for a given protein. `kind` is restricted to a
    fixed whitelist so this can't be used to read arbitrary files."""
    if kind not in _DOWNLOAD_KINDS:
        raise Http404("Unknown file kind.")

    folder_name = sanitize_folder_name(protein_name)
    filename_func, content_type = _DOWNLOAD_KINDS[kind]
    file_path = os.path.join(settings.PIPELINE_OUTPUT_DIR, folder_name, filename_func(folder_name))
    if not os.path.exists(file_path):
        raise Http404("File not found.")

    return FileResponse(open(file_path, "rb"), content_type=content_type,
                         filename=os.path.basename(file_path))
