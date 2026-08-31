import hashlib
import os
import re
import shutil
import sys
import traceback

from django.conf import settings
from django.core.paginator import Paginator
from django.http import FileResponse, Http404
from django.shortcuts import get_object_or_404, render, redirect
from django.urls import reverse
from django.utils import timezone
from django.views.decorators.http import require_POST

from .forms import SequenceSearchForm
from .models import SearchRun
from .data import (
    build_peg_plot, build_temp_plot, build_table_rows, build_method_legend,
    load_conditions, score_gradient_stops, TABLE_COLUMNS, SCORE_MIN, SCORE_MAX, SCORE_TICKS,
)

# Make the pipeline package importable (pipeline/ sits at the project root).
sys.path.insert(0, str(settings.BASE_DIR))
from pipeline.pdb_sequence_search import find_homologs_with_conditions, filter_experimental_conditions
from pipeline.compound_extraction import process_csv
from pipeline.plot import run_plot


# Real project metadata (kept in sync with CITATION.cff / pcde_metadata.yaml)
# for the Organization page -- not invented copy.
LAB_ORGANIZATION = {
    "name": "RitAreaSciencePark",
    "type": "Research Organization",
    "country": "Italy",
    "github_url": "https://github.com/RitAreaSciencePark",
}
LAB_TEAM = [
    {
        "name": "Ruth Nana Njantang",
        "orcid": "0000-0002-6003-7521",
    },
    {
        "name": "Valerio Piomponi",
        "orcid": "0000-0003-0433-8319",
    },
    {
        "name": "Andrea Dalle Vedove",
        "orcid": "0000-0001-5127-7737",
    },
]


def _normalize_sequence(sequence: str) -> str:
    """Collapse whitespace and case so repeated searches with the same
    biological sequence are treated as the same search."""
    if not sequence:
        return ""
    return "".join(sequence.strip().split()).upper()


def _build_search_signature(protein_name: str, sequence: str, sequence_type: str,
                           identity: float, evalue: float, max_hits: int,
                           llm_fallback: bool) -> str:
    payload = "|".join([
        (protein_name or "").strip(),
        _normalize_sequence(sequence),
        (sequence_type or "").strip(),
        str(float(identity)),
        str(float(evalue)),
        str(int(max_hits)),
        "1" if llm_fallback else "0",
    ])
    return hashlib.sha256(payload.encode("utf-8")).hexdigest()


def _existing_completed_run(**filters):
    return SearchRun.objects.filter(**filters).order_by("-finished_at", "-created_at").first()


def home(request):
    """Landing page: what Crystal Explorer is and why it exists."""
    return render(request, "viewer/home.html", {"active_nav": "home"})


def organization(request):
    """Who's behind the tool -- organization + team, from the repo's own
    citation/metadata files rather than separately-maintained copy."""
    return render(request, "viewer/Organization.html", {
        "active_nav": "organization",
        "organization": LAB_ORGANIZATION,
        "team": LAB_TEAM,
    })


def sanitize_folder_name(name: str) -> str:
    """Same sanitizer as the CLI's main.py, kept in sync here."""
    name = name.strip().replace(" ", "_")
    name = re.sub(r"[^A-Za-z0-9_\-]", "", name)
    return name or "protein_output"


def unique_folder_name(protein_name: str) -> str:
    """Timestamp-suffixed output folder name for one search run.

    Two searches for the same protein name used to collide on the same
    folder -- a rerun (different sequence, different thresholds, or just
    a repeat months later) silently overwrote the previous run's outputs,
    which defeats being able to keep a history of past runs. Each run now
    gets its own folder, so history entries always resolve to the exact
    files that produced them."""
    base = sanitize_folder_name(protein_name)
    stamp = timezone.now().strftime("%Y%m%d-%H%M%S")
    candidate = f"{base}_{stamp}"
    suffix = 1
    unique = candidate
    while os.path.exists(os.path.join(settings.PIPELINE_OUTPUT_DIR, unique)):
        unique = f"{candidate}-{suffix}"
        suffix += 1
    return unique


def index(request):
    """GET: show the search form. POST: run the full pipeline (save FASTA
    -> search -> filter -> compound extraction -> plots/table)
    synchronously, then redirect to the results page for that protein.

    Repeated searches with the same protein sequence and threshold settings
    are now deduplicated against SearchRun.search_signature: the existing
    history row is reused, and the app redirects straight to that result
    instead of inserting a second row.
    """
    if request.method == "POST":
        form = SequenceSearchForm(request.POST)
        if form.is_valid():
            protein_name = form.cleaned_data["protein_name"]
            sequence = form.cleaned_data["sequence"]
            normalized_sequence = _normalize_sequence(sequence)
            signature = _build_search_signature(
                protein_name,
                normalized_sequence,
                form.cleaned_data["sequence_type"],
                form.cleaned_data["identity"],
                form.cleaned_data["evalue"],
                form.cleaned_data["max_hits"],
                form.cleaned_data["llm_fallback"],
            )

            existing_run = SearchRun.objects.filter(search_signature=signature).first()
            if existing_run:
                return redirect(reverse("viewer:results", args=[existing_run.folder_name]))

            folder_name = unique_folder_name(protein_name)
            output_dir = os.path.join(settings.PIPELINE_OUTPUT_DIR, folder_name)
            os.makedirs(output_dir, exist_ok=True)
            started_at = timezone.now()

            run = SearchRun.objects.create(
                protein_name=protein_name,
                folder_name=folder_name,
                sequence=normalized_sequence,
                search_signature=signature,
                sequence_type=form.cleaned_data["sequence_type"],
                identity=form.cleaned_data["identity"],
                evalue=form.cleaned_data["evalue"],
                max_hits=form.cleaned_data["max_hits"],
                llm_fallback=form.cleaned_data["llm_fallback"],
                sequence_preview=normalized_sequence[:80],
            )

            try:
                # Save the input sequence as a FASTA file in the protein's
                # output folder -- this is also offered as a static export
                # on the results page.
                fasta_path = os.path.join(output_dir, f"{folder_name}.fasta")
                with open(fasta_path, "w") as f:
                    f.write(f">{folder_name}\n{normalized_sequence}\n")

                rows = find_homologs_with_conditions(
                    normalized_sequence,
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
                run.status = SearchRun.STATUS_FAILED
                run.error_message = str(e)
                run.finished_at = timezone.now()
                run.runtime_seconds = (run.finished_at - started_at).total_seconds()
                run.save()
                return render(request, "viewer/explorer.html", {
                    "form": form,
                    "run_error": f"Pipeline failed: {e}",
                    "run_error_detail": traceback.format_exc(),
                    "active_nav": "explorer",
                })

            run.status = SearchRun.STATUS_COMPLETED
            run.row_count = _grouped_row_count(output_dir)
            run.finished_at = timezone.now()
            run.runtime_seconds = (run.finished_at - started_at).total_seconds()
            run.save()

            return redirect(reverse("viewer:results", args=[folder_name]))
    else:
        form = SequenceSearchForm()

    return render(request, "viewer/explorer.html", {"form": form, "active_nav": "explorer"})


def _grouped_row_count(output_dir: str):
    """Row count of Grouped_conditions.csv, if it was produced -- used to
    show "N conditions found" on the history list without re-reading the
    full CSV/plots every time that page is loaded."""
    path = os.path.join(output_dir, "Grouped_conditions.csv")
    if not os.path.exists(path):
        return None
    try:
        with open(path) as f:
            return max(sum(1 for _ in f) - 1, 0)  # minus the header row
    except OSError:
        return None


def history(request):
    """Past searches, newest first, so a lab can revisit a run (or see why
    one failed) without re-submitting the form."""
    runs = SearchRun.objects.all()
    paginator = Paginator(runs, 25)
    page_number = request.GET.get("page")
    page_obj = paginator.get_page(page_number)

    return render(request, "viewer/history.html", {"page_obj": page_obj, "active_nav": "history"})


@require_POST
def delete_run(request, run_id):
    """Remove one history entry -- POST-only (the template's delete button
    confirms first) so it can't be triggered by a stray GET/crawler. Also
    removes the run's output folder, since a history row with no files
    behind it is more confusing than useful."""
    run = get_object_or_404(SearchRun, id=run_id)

    output_dir = os.path.join(settings.PIPELINE_OUTPUT_DIR, run.folder_name)
    if os.path.isdir(output_dir):
        shutil.rmtree(output_dir, ignore_errors=True)
    run.delete()

    page = request.POST.get("page")
    redirect_url = reverse("viewer:history")
    if page:
        redirect_url += f"?page={page}"
    return redirect(redirect_url)


def results(request, protein_name):
    """Show the plots + table for a protein that's already been run."""
    run_record = SearchRun.objects.filter(folder_name=protein_name).first()
    if not run_record:
        run_record = SearchRun.objects.filter(protein_name=protein_name).order_by("-finished_at", "-created_at").first()

    output_dir = os.path.join(settings.PIPELINE_OUTPUT_DIR, protein_name)
    grouped_csv_path = os.path.join(output_dir, "Grouped_conditions.csv")

    error = None
    peg_plot_div = None
    temp_plot_div = None
    table_rows = []
    method_legend = []
    run_runtime_seconds = run_record.runtime_seconds if run_record else None

    if not os.path.exists(grouped_csv_path):
        error = (
            f"No results found for '{protein_name}' yet "
            f"(expected {grouped_csv_path}). Run a search first."
        )
    else:
        try:
            import plotly.offline as pyo
            df = load_conditions(grouped_csv_path)

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
        "active_nav": "explorer",
        "form": SequenceSearchForm(initial={"protein_name": protein_name}),
        "protein_name": protein_name,
        "run_runtime_seconds": run_runtime_seconds,
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
