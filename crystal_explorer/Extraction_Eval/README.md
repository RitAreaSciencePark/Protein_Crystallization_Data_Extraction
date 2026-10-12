# Evaluation reproducibility package (ACDSA 2027 paper, Section VII)

This measures how effective `compound_extraction.py` (the dictionary +
regex extractor, no LLM involved) is at pulling structured
(compound, concentration) pairs out of real `pdbx_details` text, against
a manually curated ground truth -- precision, recall, F1, and entry-level
accuracy, the same structure as a standard extraction benchmark.

**This directory now holds the paper's current, reported result: the
665-entry held-out evaluation.** The files here are direct copies of
`heldout/heldout_*_665.json` -- see `heldout/README.md` for the full
methodology, the 100-protein panel it was drawn from, and the earlier
100-, 500-, and 1,000-entry draws it superseded. The 665-entry panel is
the *entire* remaining unused portion of that pool (not a round number
by choice -- a disjoint 1,000-entry draw is no longer possible once the
100+500+1,000-entry panels are excluded), so the pool is now fully
exhausted: every retained hit it ever produced has been manually
audited at least once across the four panels combined.

> **Note on this folder's history:** this directory originally held a
> different, smaller *development* panel (8 proteins, 60 entries,
> `seed=7`), iteratively re-examined across five rounds of dictionary
> fixes -- narrated in full below under "Five fix rounds". That panel's
> manually-judged ground truth (`ground_truth_annotations.json`, with a
> per-entry `issue` note explaining each judgment call) was overwritten
> in place when this folder was updated to show the 1,000-entry result
> instead, and it was not backed up first -- it is not recoverable. The
> raw per-protein pipeline output those 60 entries were drawn from
> (`runs/<protein>/Output_compounds.csv`, 8 proteins) is still intact
> below, and `eval_targets.json`/`eval_results.json` are trivially
> re-fetched, but the specific manual TP/FP/FN judgments themselves are
> gone -- so the fix-round narrative is kept as historical documentation
> of how the dictionary was actually developed, not as a currently
> reproducible dataset. Re-annotating all 60 entries from scratch would
> produce new judgments, not a reconstruction of the originals.

## Reproduce (current 665-entry result)

```bash
cd crystal_explorer
.venv/bin/python3 Extraction_Eval/heldout/fetch_targets.py          # -> heldout/heldout_targets.json (100 proteins)
.venv/bin/python3 Extraction_Eval/heldout/run_heldout_eval.py        # -> heldout/heldout_results.json, heldout/runs/<protein>/Output_compounds.csv (slow: live RCSB calls)
.venv/bin/python3 Extraction_Eval/heldout/build_sample_665.py        # -> heldout/heldout_sample_665_raw.json (665 entries -- every entry not already used in the 100-, 500-, or 1,000-entry draws)
# --- manual step: single-annotator audit of all 665 entries, split into
#     14 batches (thirteen of 50, one of 15; see heldout/README.md and
#     heldout/batches_665/) ---
# merge the 14 batch annotation files, then:
python3 -c "import json; d=[e for i in range(1,15) for e in json.load(open(f'Extraction_Eval/heldout/batches_665/batch_{i:02d}_annotations.json'))]; json.dump(d, open('Extraction_Eval/heldout/heldout_annotations_665.json','w'), indent=2)"
# then compute the summary (precision/recall/F1/entry-exact) from that merged file
```

The files in *this* directory (`eval_targets.json`, `eval_results.json`,
`ground_truth_sample_raw.json`, `ground_truth_annotations.json`,
`eval_summary.json`) are then just `cp`'d from the corresponding
`heldout/heldout_*_665.json` files above.

## Files

- `eval_targets.json` -- the 100 held-out proteins' query sequences
  (copy of `heldout/heldout_targets.json`).
- `eval_results.json` -- per-protein retained-hit counts from the live
  run (copy of `heldout/heldout_results.json`; 2,265 retained hits total
  across the 100 proteins).
- `runs/<protein>/Output_compounds.csv` -- **the original 8-protein
  development panel's** raw pipeline output (184 retained hits); still
  intact, but no longer what `ground_truth_sample_raw.json` below is
  drawn from. The 100-protein held-out panel's equivalent output lives
  at `heldout/runs/<protein>/Output_compounds.csv`.
- `ground_truth_sample_raw.json` -- the 665 entries (every entry not
  already used in the earlier 100-, 500-, or 1,000-entry draws),
  unannotated (copy of `heldout/heldout_sample_665_raw.json`).
- `ground_truth_annotations.json` -- the same 665 entries with
  manually determined ground truth and per-entry TP/FP/FN counts (copy
  of `heldout/heldout_annotations_665.json`).
- `eval_summary.json` -- the aggregate precision/recall/F1/accuracy
  numbers quoted in the paper (Table I): 99.2% precision, 93.6% recall,
  F1 0.963, 84.7% entry-level exact-match accuracy.

## Methodology

For each sampled `pdbx_details` string, the raw text is read first and
the true list of crystallization reagents it names (the "ground truth")
is determined independently of what the extractor produced. Each
ground-truth compound is then checked against
`compound_extraction.py`'s actual output for that entry:

- **True positive (TP):** a ground-truth compound the extractor
  correctly identified (by canonical name/chemical family).
- **False negative (FN):** a ground-truth compound the extractor missed
  entirely.
- **False positive (FP):** a compound the extractor reported that does
  not correspond to any real reagent in the text (e.g. a wrong,
  over-generalized canonical name).

Precision = TP/(TP+FP), Recall = TP/(TP+FN), F1 = 2PR/(P+R). A separate,
stricter **entry-level exact-match accuracy** counts an entry as correct
only if every compound in it was extracted with no misses and no false
positives (FN=0 and FP=0 for that entry).

Scoring is at the entity level, not the literal-mention level: if the
same true compound is named more than once in one entry (e.g. a full name
followed by its parenthetical abbreviation, ``bromophenol blue (BPB)'',
which `match_clause_all()` can match twice within one clause), that
counts as a single TP, not one TP plus spurious FPs for the repeats. This
mirrors how comparable extraction benchmarks (including aTUNApy's) score
per distinct entity per document, not per string occurrence.

## Limitations of this evaluation

This is a single-annotator manual audit, not an independently
cross-validated or crowd-sourced ground truth -- stated as such in the
paper. Borderline judgment calls (e.g. whether a bare ionic species like
"Zn(II)" or an unqualified "peg" with no molecular weight counts as a
"ground truth compound" the extractor should have resolved) are
documented per-entry in `ground_truth_notes` rather than silently
decided.

## Five fix rounds (historical -- development panel, no longer backed by data in this folder)

The following narrates how `compound_extraction.py`'s dictionary was
actually built: five rounds of fixes driven by the original 8-protein,
60-entry development panel, re-examining the same panel after each
change. **The manually-judged data this narrative was computed from no
longer exists in this repository** (see the note at the top of this
file) -- the numbers and entry IDs below are kept as a historical record
of the tool's development, not as something you can currently
reproduce by re-running a script here.

1. **Implementation bugs** (found before the ground-truth audit, fixed
   first): `_format_concentration()` silently dropped correctly-resolved
   non-PEG concentration ranges to `"nan"`; `match_clause()`'s
   one-match-per-clause design silently dropped a second reagent from
   clauses `split_into_clauses()` failed to separate. Both fixed via a
   new `match_clause_all()` and a range-aware `_format_concentration()`.
2. **Missing alias spellings** (found *by* the ground-truth audit): the
   first 60-entry run scored 130 TP / 5 FP / 15 FN (precision 0.963,
   recall 0.897, F1 0.929); three of those errors (entries `3ED9`,
   `2ZSR`, `2BLP`) traced to two missing alias spellings --
   `"trishcl"` (no separator) for Tris-HCl, and the reversed word order
   `"peg 5000 mme"` for PEG MME 5000. Both added as one-line dictionary
   entries. Re-running the same evaluation confirmed the fix: 133 TP /
   4 FP / 12 FN (precision 0.971, recall 0.917, F1 0.943).
3. **Tail-token mislabeling**: the 4 FP / 4 FN left after round 2
   (entries `1IF7`, `1IF8`) were both the same pattern -- ``methyl
   mercuric acetate`` and ``tris-sulfate`` collapsing to the generic
   ``acetate``/``tris`` aliases. Registering both as their own canonical
   dictionary entries (`REAGENT_DICTIONARY["methyl mercuric acetate"]`
   and `REAGENT_DICTIONARY["tris-sulfate"]`, plus matching
   `REAGENT_CATEGORY` entries) resolved both, since aliases are matched
   longest-first. Re-running the same evaluation confirmed it: 137 TP /
   0 FP / 8 FN (precision 1.000, recall 0.945, F1 0.972) -- zero false
   positives from this point on.
4. **Remaining coverage gaps** (found *by* the ground-truth audit, the
   8 FN left after round 3): added `"magnesium citrate"`, `"PEG 1550"`,
   `"zn(ch2h3o2)2"` as a zinc-acetate alias, `"Bromophenol blue"`
   (aliases `"bromophenol blue"`, `"bromo phenyl blue"` -- the exact
   typo'd variant in the source text -- and `"bpb"`), and a generic
   `"zinc"` entry (`"zn(ii)"` alias) for a bare, unspecified counter-ion
   mention. Also added `"ches buffer"` as a CHES alias (fixes a
   concentration-attachment miss, not an FN) and
   `"na-acetate"` as a sodium-acetate alias (sharpens two
   already-correct-but-generalized matches, not an FN either). Six of
   the eight FN closed this round; the remaining two (`"AS"` for
   ammonium sulfate, bare `"peg"` with no molecular weight) were
   deliberately left unresolved, since a dictionary cannot safely add
   such short, generic tokens as aliases without risking false matches
   on unrelated text elsewhere. Re-running the same evaluation a fourth
   time confirmed it: 143 TP / 0 FP / 2 FN (precision 1.000, recall
   0.986, F1 0.993).

   This round also fixed a pre-existing data bug found incidentally
   while adding the zinc-acetate alias: `REAGENT_DICTIONARY` defined the
   key `"zinc acetate"` twice (once with 3 aliases, once with just 1);
   Python dict-literal semantics let the second definition silently
   overwrite the first, discarding the `"zn(oac)2"` alias. The two were
   merged into one entry.

   The zinc fix also surfaced a ground-truth correction, unrelated to
   the code: entry `1C2G` actually names *two* distinct Zn(II) soak
   steps at different concentrations (1.0mM and 46mM), not one -- the
   original manual annotation had only recorded the first. Once the
   extractor started recovering both, the ground truth was corrected to
   match the source text: **144 TP / 0 FP / 2 FN (precision 1.000,
   recall 0.986, F1 0.993)** after this correction.
5. **Parsing and unit gaps** (found *by* the ground-truth audit, the
   2 FN left after round 4): added `"mmolar"`/`"millimolar"`/`"molar"`
   as recognized concentration units (normalized to mM/M), recovering a
   previously-unrecognized `"200 MMOLAR"` concentration. Changed
   `split_into_clauses()` to normalize embedded line-wrap newlines to
   spaces instead of splitting on them -- PDB deposition text is
   typically fixed-width wrapped mid-sentence in the mmCIF file, not
   genuinely multi-clause at every line break -- recovering a
   magnesium-sulfate concentration a stray newline had severed from its
   amount. Added `"as"` as an ammonium-sulfate alias, but restricted it
   (via a new `_CONTEXT_ONLY_ALIASES` set) to only match in patterns
   A/B, immediately adjacent to a stated amount and unit, never as a
   bare pattern-C mention -- the word "as" is far too common in ordinary
   English narrative to trust standalone. Making that fix take effect on
   its actual target case (`2EMN`, `"WELL SOLUTION: 2.1 M AS"`) also
   required loosening the pH-clause filter, which had been discarding
   the whole clause because it happened to start with a pH restatement
   even though real reagent data followed. Only the `"as"` fix closed a
   false negative (1 of the 2 remaining); the unit and newline fixes
   each recovered a previously-lost concentration on an
   already-correctly-identified compound. Re-running the same
   evaluation a fifth time, and separately checking every other sampled
   entry for regressions, confirmed it: **145 TP / 0 FP / 1 FN
   (precision 1.000, recall 0.993, F1 0.997)** -- this panel's final
   numbers. The one remaining miss (bare `"peg"`, no molecular weight)
   was left as a deliberate non-fix.

The dictionary that resulted from these five rounds is what the
1,000-entry held-out evaluation was originally run against. That panel
was then, in effect, used as a sixth development round -- an extended
session of dictionary and parsing fixes driven by manually reading
entries from it -- after which the 665-entry panel above was run as the
genuinely fresh held-out check (see `heldout/README.md`), with no
further tuning based on what the 665-entry panel found.
