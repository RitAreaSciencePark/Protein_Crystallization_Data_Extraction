# Crystal Explorer

A Django app that runs the full crystallization-condition pipeline
straight from a web form:

1. **You type/paste** a protein name and a sequence (raw or FASTA) into
   the form.
2. The app runs, server-side, in one request:
   - `pipeline/pdb_sequence_search.py` -- searches RCSB PDB (X-ray only)
     for sequence homologs and pulls crystallization metadata.
   - `pipeline/compound_extraction.py` -- extracts compounds/concentrations
     from `pdbx_details`.
   - `pipeline/plot.py` -- builds the pH-vs-Temperature and pH-vs-PEG
     plots, a colored PDF summary table, and `Grouped_conditions.csv`.
3. You're redirected to a results page: an **interactive** Plotly version
   of the pH-vs-PEG plot next to the full table -- click a point and it
   scrolls to and highlights the matching table row. Static PDF/PNG
   versions are available to download alongside it.

Each protein's run gets its own folder under `pipeline_outputs/<name>/`,
so you can search a new protein or come back to `/protein/<name>/` for
one you already ran.

## Setup

```bash
pip install -r requirements.txt
python manage.py migrate
python manage.py runserver
```

Then open http://127.0.0.1:8000/, fill in the form, and hit **Run search**.

To use the LLM fallback for compound extraction (resolves reagents the
built-in dictionary doesn't recognize, and learns them for next time),
set `ANTHROPIC_API_KEY` in the environment before starting the server,
and check the "Use Claude..." box in the form.

## Notes on running the search synchronously

The search view runs stage 1-3 of the pipeline **inside the HTTP
request** -- simplest to set up, but it means the browser tab waits
for the whole thing (RCSB search + per-hit mmCIF fetches + plotting) to
finish before the results page loads. The form defaults `max_hits` to 25
to keep this reasonable. For heavier use (many hits, many users), moving
this to a background task (Celery, Django-Q, etc.) that the results page
polls for would be the next step -- ask if you want that built out.

## Project layout

```
config/                  Django settings/urls
pipeline/                 the 3-stage pipeline, copied in as an importable package
  pdb_sequence_search.py
  compound_extraction.py
  plot.py
viewer/
  forms.py                protein name + sequence input form
  views.py                index (form/run) + results + download views
  data.py                 Grouped_conditions.csv -> Plotly figure + table rows
  templates/viewer/
    base.html
    index.html             the input form
    results.html            plot + table + click-to-highlight JS
pipeline_outputs/          created at runtime, one subfolder per protein
```
