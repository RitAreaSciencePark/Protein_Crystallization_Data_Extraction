import hashlib
import os
import re
import shutil
import sys
import threading
import traceback

import requests

from django.conf import settings
from django.contrib import messages
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
    build_reagent_frequency, build_reagent_frequency_plot, build_reagent_frequency_by_category_plot,
    build_summary_stats,
    load_conditions, score_gradient_stops, TABLE_COLUMNS, SCORE_MIN, SCORE_MAX, SCORE_TICKS,
)

# Make the pipeline package importable (pipeline/ sits at the project root).
sys.path.insert(0, str(settings.BASE_DIR))
from pipeline.pdb_sequence_search import find_homologs_with_conditions, filter_experimental_conditions
from pipeline.compound_extraction import process_csv
from pipeline.plot import run_plot
from pipeline.metadata_generator import generate_metadata_for_outputs


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
                           identity: float, evalue: float, max_hits: int) -> str:
    payload = "|".join([
        (protein_name or "").strip().lower(),
        _normalize_sequence(sequence),
        (sequence_type or "").strip(),
        str(float(identity)),
        str(float(evalue)),
        str(int(max_hits)),
    ])
    return hashlib.sha256(payload.encode("utf-8")).hexdigest()


def _build_sequence_signature(protein_name: str, sequence: str, sequence_type: str) -> str:
    """Signature of just the search *subject* -- protein name, sequence,
    and sequence type -- deliberately excluding the identity/E-value/
    max-hits thresholds and the LLM-fallback flag. Two searches sharing
    this signature are "the same search with different parameters": the
    later one reuses and overwrites the earlier run's history row and
    output folder instead of getting its own (see index()).

    The protein name is lowercased before hashing (matching
    _build_search_signature) so "PLpro" and "plpro" are treated as the
    same search -- only the sequence itself has to be identical."""
    payload = "|".join([
        (protein_name or "").strip().lower(),
        _normalize_sequence(sequence),
        (sequence_type or "").strip(),
    ])
    return hashlib.sha256(payload.encode("utf-8")).hexdigest()


def _friendly_pipeline_error(e: Exception) -> str:
    """RCSB's search API answers a malformed query (e.g. a sequence with
    characters that aren't valid for the chosen protein/DNA/RNA type, or
    one that's empty/too short) with a 400 -- surfaced by `requests` as an
    HTTPError whose message is just "400 Client Error: Bad Request for
    url: ...". That's meaningless to someone who just typed a sequence in,
    so translate it into something actionable instead of the raw
    exception text."""
    if isinstance(e, requests.exceptions.HTTPError) and e.response is not None \
            and e.response.status_code == 400:
        return (
            "Bad format for sequence -- please check that it only contains "
            "valid characters for the selected sequence type (protein, DNA, "
            "or RNA) and try again."
        )
    return f"Pipeline failed: {e}"


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


def _current_session_key(request) -> str:
    """The visitor's anonymous Django session key, creating the session
    (and its cookie) on the spot if this is their first request -- Django
    doesn't allocate a key until something is stored in the session. No
    login, no IP: just the same per-browser cookie Django already sets."""
    if not request.session.session_key:
        request.session.save()
    return request.session.session_key


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


def _set_progress(run: "SearchRun", percent: int, label: str) -> None:
    """Update the run's coarse progress both on the in-memory `run` (so the
    final run.save() in _run_pipeline_in_background doesn't clobber it back
    to 0) and directly in the DB via a targeted .update() (so the results()
    view polling from another request/thread sees it immediately without
    waiting for that final save)."""
    run.progress_percent = percent
    run.progress_label = label
    SearchRun.objects.filter(id=run.id).update(progress_percent=percent, progress_label=label)


def _run_pipeline_in_background(run_id: int, output_dir: str, folder_name: str,
                                 normalized_sequence: str, sequence_type: str,
                                 identity: float, evalue: float, max_hits: int,
                                 lab_name: str, started_at) -> None:
    """The actual search->filter->extract->plot->metadata pipeline, run on a
    background thread so the request that triggered it can redirect
    immediately instead of blocking for the whole run (RCSB search +
    extraction can take anywhere from seconds to a couple of minutes).

    Runs in its own OS thread, not the request thread, so it re-fetches the
    SearchRun by id (an ORM instance from the request thread isn't safe to
    share across threads) and closes its own DB connections when done --
    standard practice for Django code running outside the request/response
    cycle, so connections don't pile up across many searches.

    This is deliberately a plain background thread, not a real task queue
    (Celery/django-q/huey): there's no separate worker process to deploy or
    keep alive, which matches how this app is actually run today (`manage.py
    runserver`, nothing else). The real tradeoff against a proper queue: if
    the server process restarts mid-search, this thread dies with it and the
    run is left stuck at "running" with no automatic retry -- acceptable for
    a single-process lab tool, not for a multi-worker production deployment,
    which is exactly the "natural next step" the paper already flags this as
    one layer short of."""
    from django.db import close_old_connections

    close_old_connections()
    run = SearchRun.objects.get(id=run_id)
    try:
        _set_progress(run, 5, "Preparing input sequence")
        fasta_path = os.path.join(output_dir, f"{folder_name}.fasta")
        with open(fasta_path, "w") as f:
            f.write(f">{folder_name}\n{normalized_sequence}\n")

        _set_progress(run, 10, "Searching RCSB PDB for homologs")
        rows = find_homologs_with_conditions(
            normalized_sequence,
            identity_cutoff=identity,
            evalue_cutoff=evalue,
            sequence_type=sequence_type,
            max_hits=max_hits,
            verbose=True,
        )

        _set_progress(run, 60, "Filtering experimental conditions")
        out_path = os.path.join(output_dir, "Output.csv")
        filter_experimental_conditions(rows, output_csv=out_path)

        _set_progress(run, 70, "Extracting compounds from crystallization conditions")
        compounds_out_path = os.path.join(output_dir, "Output_compounds.csv")
        process_csv(out_path, compounds_out_path)

        _set_progress(run, 85, "Building plots and summary table")
        run_plot(output_dir, folder_name)

        _set_progress(run, 95, "Writing FAIR metadata")
        generate_metadata_for_outputs(
            folder_name, output_dir,
            sequence_type=sequence_type,
            search_params={
                "identity": identity,
                "evalue": evalue,
                "max_hits": max_hits,
                "lab_name": lab_name,
            },
        )
    except Exception as e:
        # Printed server-side (visible in the runserver/WSGI log) since
        # there's no request left to render a detailed traceback into --
        # the friendly, one-line version is what a user sees on the
        # results page via run.error_message.
        traceback.print_exc()
        run.status = SearchRun.STATUS_FAILED
        run.error_message = _friendly_pipeline_error(e)
        run.finished_at = timezone.now()
        run.runtime_seconds = (run.finished_at - started_at).total_seconds()
        run.save()
        close_old_connections()
        return

    run.status = SearchRun.STATUS_COMPLETED
    run.progress_percent = 100
    run.progress_label = ""
    run.row_count = _grouped_row_count(output_dir)
    run.finished_at = timezone.now()
    run.runtime_seconds = (run.finished_at - started_at).total_seconds()
    run.save()
    close_old_connections()


def index(request):
    """GET: show the search form. POST: kick off the full pipeline (save
    FASTA -> search -> filter -> compound extraction -> plots/table) on a
    background thread and redirect immediately to the results page, which
    shows a "still running" state until the thread finishes -- see
    _run_pipeline_in_background() for why a thread and not a real queue.

    Repeated searches are deduplicated at two levels:
      - Same protein/sequence/type AND same thresholds/flags (full
        search_signature match): the existing completed run is reused
        outright -- redirect straight to it, no re-run.
      - Same protein/sequence/type but different thresholds/flags
        (sequence_signature match only): treated as "redo this search
        with different parameters" -- the existing history row and output
        folder are reused and overwritten in place, rather than piling up
        a separate history entry per parameter tweak.
      - No match at all: a brand-new history row and output folder.
    """
    if request.method == "POST":
        form = SequenceSearchForm(request.POST)
        if form.is_valid():
            protein_name = form.cleaned_data["protein_name"]
            sequence = form.cleaned_data["sequence"]
            normalized_sequence = _normalize_sequence(sequence)
            sequence_type = form.cleaned_data["sequence_type"]
            signature = _build_search_signature(
                protein_name,
                normalized_sequence,
                sequence_type,
                form.cleaned_data["identity"],
                form.cleaned_data["evalue"],
                form.cleaned_data["max_hits"],
            )

            existing_run = SearchRun.objects.filter(search_signature=signature).first()
            if existing_run:
                return redirect(reverse("viewer:results", args=[existing_run.folder_name]))

            sequence_signature = _build_sequence_signature(protein_name, normalized_sequence, sequence_type)
            reusable_run = SearchRun.objects.filter(sequence_signature=sequence_signature).first()

            if reusable_run:
                run = reusable_run
                folder_name = run.folder_name
                output_dir = os.path.join(settings.PIPELINE_OUTPUT_DIR, folder_name)
                os.makedirs(output_dir, exist_ok=True)
                run.search_signature = signature
                run.status = SearchRun.STATUS_RUNNING
                run.progress_percent = 0
                run.progress_label = ""
            else:
                folder_name = unique_folder_name(protein_name)
                output_dir = os.path.join(settings.PIPELINE_OUTPUT_DIR, folder_name)
                os.makedirs(output_dir, exist_ok=True)
                run = SearchRun(folder_name=folder_name, search_signature=signature)

            started_at = timezone.now()
            run.protein_name = protein_name
            run.session_key = _current_session_key(request)
            run.sequence = normalized_sequence
            run.sequence_signature = sequence_signature
            run.sequence_type = sequence_type
            run.identity = form.cleaned_data["identity"]
            run.evalue = form.cleaned_data["evalue"]
            run.max_hits = form.cleaned_data["max_hits"]
            run.lab_name = form.cleaned_data["lab_name"]
            run.sequence_preview = normalized_sequence[:80]
            run.error_message = ""
            run.save()

            # Remember the lab name in this browser's session so the next
            # search pre-fills it -- a convenience, not an account; an
            # empty value here just means the field stays blank next time.
            if run.lab_name:
                request.session["lab_name"] = run.lab_name

            thread = threading.Thread(
                target=_run_pipeline_in_background,
                args=(
                    run.id, output_dir, folder_name, normalized_sequence, sequence_type,
                    form.cleaned_data["identity"], form.cleaned_data["evalue"],
                    form.cleaned_data["max_hits"], form.cleaned_data["lab_name"], started_at,
                ),
                daemon=True,
            )
            thread.start()

            return redirect(reverse("viewer:results", args=[folder_name]))
    else:
        remembered_lab = request.session.get("lab_name", "")
        form = SequenceSearchForm(initial={"lab_name": remembered_lab} if remembered_lab else None)

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
    """This visitor's own past searches, newest first, so they can revisit
    a run (or see why one failed) without re-submitting the form. Scoped
    to the current anonymous session -- not an account, just the same
    cookie used for the delete-ownership check below. This is independent
    of the dedup in index(), which still checks every SearchRun regardless
    of session, so a repeated search still reuses the existing run instead
    of re-running the pipeline even though this list only shows one
    visitor's slice of that shared data."""
    runs = SearchRun.objects.filter(session_key=_current_session_key(request))
    paginator = Paginator(runs, 25)
    page_number = request.GET.get("page")
    page_obj = paginator.get_page(page_number)

    return render(request, "viewer/history.html", {
        "page_obj": page_obj,
        "active_nav": "history",
    })


@require_POST
def delete_run(request, run_id):
    """Remove one history entry -- POST-only (the template's delete button
    confirms first) so it can't be triggered by a stray GET/crawler. Also
    removes the run's output folder, since a history row with no files
    behind it is more confusing than useful.

    Only the browser session that (most recently) ran a search may delete
    it -- the history list itself stays visible to everyone, but deleting
    someone else's row is blocked. A row with no session_key at all (from
    before this check existed) has no owner to enforce, so it stays
    deletable by anyone rather than becoming permanently stuck."""
    run = get_object_or_404(SearchRun, id=run_id)

    page = request.POST.get("page")
    redirect_url = reverse("viewer:history")
    if page:
        redirect_url += f"?page={page}"

    if run.session_key and run.session_key != _current_session_key(request):
        messages.error(
            request,
            f"“{run.protein_name}” was submitted from a different "
            "browser session, so it can only be deleted from the browser "
            "that ran it.",
        )
        return redirect(redirect_url)

    output_dir = os.path.join(settings.PIPELINE_OUTPUT_DIR, run.folder_name)
    if os.path.isdir(output_dir):
        shutil.rmtree(output_dir, ignore_errors=True)
    run.delete()

    return redirect(redirect_url)


def results(request, protein_name):
    """Show the plots + table for a protein that's already been run."""
    run_record = SearchRun.objects.filter(folder_name=protein_name).first()
    if not run_record:
        run_record = SearchRun.objects.filter(protein_name=protein_name).order_by("-finished_at", "-created_at").first()

    output_dir = os.path.join(settings.PIPELINE_OUTPUT_DIR, protein_name)
    grouped_csv_path = os.path.join(output_dir, "Grouped_conditions.csv")

    error = None
    no_hits = False
    running = False
    progress_percent = 0
    progress_label = ""
    peg_plot_div = None
    temp_plot_div = None
    reagent_plot_div = None
    reagent_category_plot_div = None
    table_rows = []
    method_legend = []
    summary_stats = None
    run_runtime_seconds = run_record.runtime_seconds if run_record else None

    if run_record and run_record.status == SearchRun.STATUS_RUNNING:
        # The pipeline is still running on its background thread (see
        # _run_pipeline_in_background) -- don't even look for the CSV yet,
        # since a partially-written file mid-run is worse than none. The
        # template auto-refreshes this page until the status changes.
        running = True
        progress_percent = run_record.progress_percent
        progress_label = run_record.progress_label
    elif run_record and run_record.status == SearchRun.STATUS_FAILED:
        error = run_record.error_message or f"Search for '{protein_name}' failed."
    elif not os.path.exists(grouped_csv_path):
        error = (
            f"No results found for '{protein_name}' yet "
            f"(expected {grouped_csv_path}). Run a search first."
        )
    else:
        try:
            import plotly.offline as pyo
            df = load_conditions(grouped_csv_path)

            if df.empty:
                # The search ran fine but matched zero PDB entries (or none
                # of the matches had usable crystallization conditions) --
                # that's a normal outcome, not an error, so show a plain
                # "no hits" message instead of empty plots/table.
                no_hits = True
            else:
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

                reagent_freq = build_reagent_frequency(df)
                reagent_plot_div = pyo.plot(
                    build_reagent_frequency_plot(reagent_freq),
                    output_type="div", include_plotlyjs=False, config=plot_config
                )
                reagent_category_plot_div = pyo.plot(
                    build_reagent_frequency_by_category_plot(reagent_freq),
                    output_type="div", include_plotlyjs=False, config=plot_config
                )

                table_rows = build_table_rows(df)
                summary_stats = build_summary_stats(df, reagent_freq)
        except Exception as e:
            error = f"Could not load or process {grouped_csv_path}: {e}"

    pdf_path = os.path.join(output_dir, f"{protein_name}_Cryst_cocktail_Table.pdf")
    peg_png_path = os.path.join(output_dir, f"{protein_name}_PEG.png")
    temp_png_path = os.path.join(output_dir, f"{protein_name}_TEMP.png")
    fasta_path = os.path.join(output_dir, f"{protein_name}.fasta")
    compounds_csv_path = os.path.join(output_dir, "Output_compounds.csv")
    metadata_json_path = os.path.join(output_dir, f"{protein_name}_metadata.json")

    return render(request, "viewer/explorer.html", {
        "active_nav": "explorer",
        "form": SequenceSearchForm(initial={"protein_name": protein_name}),
        "protein_name": protein_name,
        "run_runtime_seconds": run_runtime_seconds,
        "peg_plot_div": peg_plot_div,
        "temp_plot_div": temp_plot_div,
        "reagent_plot_div": reagent_plot_div,
        "reagent_category_plot_div": reagent_category_plot_div,
        "method_legend": method_legend,
        "score_gradient_stops": score_gradient_stops(),
        "score_min": SCORE_MIN,
        "score_max": SCORE_MAX,
        "score_ticks": SCORE_TICKS,
        "table_rows": table_rows,
        "summary_stats": summary_stats,
        "table_columns": TABLE_COLUMNS,
        "error": error,
        "no_hits": no_hits,
        "running": running,
        "progress_percent": progress_percent,
        "progress_label": progress_label,
        "csv_path": grouped_csv_path,
        "has_pdf": os.path.exists(pdf_path),
        "has_peg_png": os.path.exists(peg_png_path),
        "has_temp_png": os.path.exists(temp_png_path),
        "has_fasta": os.path.exists(fasta_path),
        "has_compounds_csv": os.path.exists(compounds_csv_path),
        "has_metadata": os.path.exists(metadata_json_path),
    })


_DOWNLOAD_KINDS = {
    "pdf": (lambda name: f"{name}_Cryst_cocktail_Table.pdf", "application/pdf"),
    "peg": (lambda name: f"{name}_PEG.png", "image/png"),
    "temp": (lambda name: f"{name}_TEMP.png", "image/png"),
    "metadata": (lambda name: f"{name}_metadata.json", "application/json"),
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
