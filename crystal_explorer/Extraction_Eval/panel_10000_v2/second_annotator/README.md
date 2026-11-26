# Blinded second annotation of the held-out audit

Responds to the reviewer comment that the ground truth was made by a single
annotator involved in building the dictionary. The original audit was also
**not blinded**: each audit batch showed the extractor's output next to the
raw text (`../batches/batch_NN_raw.json`), which risks anchoring.

A second annotator (B) who was not involved in building the extractor
annotates a random 200 of the 1,000 audited entries (seed 2032). B sees only
the raw `pdbx_details` text.

On these 200 entries the original audit (annotator A) gives TP 573, FP 4,
FN 33: precision 0.993, recall 0.946, close to the full 1,000.

## Steps, in this order

1. **Package (done):** `python build_package.py`.
2. **Annotator B:** send B the folder `for_annotator_B/` **only**
   (`annotation_sheet.csv` + `GUIDELINES.md`). Do not send anything else
   from this directory. B returns `annotation_sheet.csv` with `reagents` filled
   in, and you put it back at `for_annotator_B/annotation_sheet.csv`.
3. **Annotator A, before seeing B's file:** review
   `annotator_A_structured.csv`. 172 rows were filled automatically from
   error-free entries (`needs_review = no`). In the 28 rows with
   `needs_review = yes`:
   - remove the extractor items you scored as false positives;
   - replace each `... [FN?]` draft with the proper reagent name, or delete
     it if it was not counted (e.g. `compound 4`);
   - add missed reagents that the notes describe without quotes (e.g. 1FXO:
     TMP).

   Then don't change the file again. The scoring script refuses to run while
   any `[FN?]` is left.
4. `python score_agreement.py prepare` writes `adjudication.csv`, listing
   every item not matched automatically.
5. **Adjudication**, ideally by A and B together: fill in `same_as` and
   `consensus` (see the docstring of `score_agreement.py`).
6. `python score_agreement.py score` writes `agreement_summary.json` and
   `disagreements.csv`.

Use the project venv (`crystal_explorer/.venv/bin/python3`) for steps 4-6,
because they need numpy.

## What gets reported

- **A vs B:** positive specific agreement (mention-level F1 between the two
  annotators) and entry-level exact agreement. Cohen's kappa is not
  defined here, since there are no true negatives.
- **Extractor scores** (precision, recall, F1, entry exact) against A
  (original), against B alone, and against the adjudicated consensus.
- **Disagreements,** listed with the adjudication decision, so they can be
  grouped by type in the paper.

All CIs come from a bootstrap over entries (2,000 resamples).

The scoring was tested end to end on a simulated annotator B in a scratch
copy. No real annotator-B data exists yet.
