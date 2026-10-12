# Held-out validation panel (ACDSA 2027 paper, Section VII-E/F)

The `../` directory (the original 8-protein / 60-entry panel) was
iteratively re-examined across five rounds of fixes to
`compound_extraction.py` -- exactly the discipline a rule-based
extractor benefits from, but also exactly the setup where a tool can
end up looking better than it really is, since the dictionary was
shaped by repeatedly looking at that same sample.

This directory is a second, **independent** evaluation: 100 proteins
with zero overlap with the original panel, sampled and manually audited
**once**, against the already-finished, five-rounds-fixed code, with
**no further code changes made in response to what it found**. It is
the more honest measurement of real-world performance.

## Result

| | Development panel (`../`) | Held-out panel (here) |
|---|---|---|
| Sampled entries | 60 | 60 |
| Precision | 100% | **97.1%** |
| Recall | 99.3% | **90.4%** |
| F1 | 0.997 | **0.936** |
| Entry-exact | 98.3% (59/60) | **80.0% (48/60)** |

The drop is expected and is the entire point of running this panel. The
4 false positives all trace to the *same* tail-token mislabeling
pattern already identified on the development panel (Section VII-C of
the paper) -- but on **new** compound pairs (`Tris-acetate`, `potassium
succinate`) neither of the two development-panel fixes anticipated,
confirming it's a real, recurring architectural weakness rather than
something fully closed by those two patches. The 14 false negatives are
mostly depositor typos and genuinely rare reagents (a membrane-protein
detergent, a spelled-out buffer name, acetic acid, phenol) -- an
open-ended long tail no fixed dictionary can fully enumerate, which is
the strongest evidence in the paper for actually wiring up the
LLM-assisted fallback.

**No code was changed in response to this panel's findings.** That's
deliberate: fixing what it found would just turn this into a sixth
development round, defeating the purpose. The findings are reported in
the paper as future work (Section VIII) instead.

## Reproduce

```bash
cd crystal_explorer
.venv/bin/python3 Extraction_Eval/heldout/fetch_targets.py         # -> heldout_targets.json (100 of ~118 candidates)
.venv/bin/python3 Extraction_Eval/heldout/run_heldout_eval.py       # -> heldout_results.json, runs/<protein>/Output_compounds.csv (slow: ~19 min, live RCSB calls)
.venv/bin/python3 Extraction_Eval/heldout/build_heldout_sample.py   # -> heldout_sample_raw.json (60 random entries, seed=11)
# --- manual step: read each of the 60 pdbx_details entries and judge the
#     extractor's output against the source text, producing
#     heldout_annotations.json (same TP/FP/FN methodology as the
#     development panel -- see ../README.md's "Methodology" section) ---
.venv/bin/python3 Extraction_Eval/heldout/analyze_heldout.py        # -> heldout_summary.json
```

## Files

- `candidates.py` -- 118 candidate (name, representative PDB ID) pairs,
  hand-selected for structural/functional diversity and explicitly
  disjoint from the 8-protein development panel. Over-provisioned since
  a handful of remembered PDB IDs turned out stale/wrong (1 of 101 tried
  before reaching 100 successes: `MAP_kinase_ERK2`/`1ERK`, HTTP 404).
- `heldout_targets.json` / `heldout_failed.json` -- the 100 resolved
  targets (deposited sequence, fetched live from RCSB's FASTA endpoint,
  not typed from memory) and the ones that failed to resolve.
- `heldout_results.json` -- per-protein retained-hit counts from the
  live run (2,265 total retained hits across all 100 proteins).
- `runs/<protein>/Output_compounds.csv` -- the actual pipeline output
  for every retained hit.
- `heldout_sample_raw.json` -- 60 entries (seed=11) randomly sampled
  from the 2,265-hit pool, unannotated.
- `heldout_annotations.json` -- the same 60 entries with manually
  determined ground truth and per-entry TP/FP/FN counts, plus a short
  note explaining every error.
- `heldout_summary.json` -- the aggregate numbers in the table above,
  computed from the file above.

## 100-entry extension (superseded by the 500-entry panel below)

The panel above was extended from 60 to 100 entries: 40 further entries
(seed=23, disjoint from the original 60) drawn from the same 2,265-hit
pool, extracted fresh against the current (already-fixed) code, and
manually audited once with the same TP/FP/FN methodology --
**no further code changes made in response to what they found**, an
even stronger held-out claim than the original 60's.

- `heldout_sample_extra40_raw.json` -- the 40 new entries (seed=23),
  unannotated, with the current extractor's fresh output attached.
- `heldout_annotations_extra40.json` -- the same 40 entries with manual
  ground truth and per-entry TP/FP/FN counts.
- `heldout_summary_100.json` -- the combined 100-entry result (97.6%
  precision, 92.8% recall, F1 0.952, 84/100 entry-exact), reported in
  an earlier draft of the paper's Table I and now superseded by the
  500-entry result below.

## 500-entry panel (superseded by the 1,000-entry panel below)

A fresh, independent n=500 sample was drawn from the same 100-protein
held-out pool (2,265 retained hits total), **explicitly excluding** the
100 entries already used above (`seed=500`, 2,165-entry remaining pool),
so this is a disjoint draw, not a superset. All 500 entries were
manually audited with the same single-annotator, read-the-raw-text-first
TP/FP/FN methodology described above, against the same already-fixed
code -- no code changes made in response to what this panel found.

Because 500 entries is far more than one sitting of manual annotation,
the audit was split into 10 batches of 50 and each batch was annotated
independently (same methodology, same instructions, same worked
calibration examples) before being merged; `batches_500/` holds every
batch's raw sample and annotations for inspection.

- `heldout_sample_500_raw.json` -- the 500 entries (seed=500),
  unannotated, with the extractor's output attached at the time.
- `heldout_annotations_500.json` -- the same 500 entries with manual
  ground truth and per-entry TP/FP/FN counts (merged from `batches_500/`).
- `heldout_summary_500.json` -- the aggregate numbers for this panel
  (97.0% precision, 87.4% recall, F1 0.919, 69.6% entry-exact), reported
  in an earlier draft of the paper's Table I and now superseded by the
  1,000-entry result below.
- `batches_500/` -- the 10 intermediate 50-entry batches
  (`batch_NN_raw.json` / `batch_NN_annotations.json`) the 500-entry
  audit was split into.
- `build_sample_500.py` -- draws the 500-entry sample, excluding every
  entry already used in the 100-entry panel.

## 1,000-entry panel (superseded by the 665-entry panel below)

A second, larger fresh sample was drawn from the same 100-protein
held-out pool, `seed=1000`, **explicitly excluding every entry already
used** in both the 100-entry panel and the 500-entry panel above
(2,165-entry remaining pool after those exclusions), so this is again a
disjoint draw. All 1,000 entries were manually audited with the same
single-annotator, read-the-raw-text-first TP/FP/FN methodology, against
the same already-fixed code -- no code changes made in response to what
this panel found. The audit was split into 20 batches of 50 entries,
each annotated independently with the same instructions and calibration
examples; `batches_1000/` holds every batch's raw sample and
annotations.

| | 100-entry panel | 500-entry panel | 1,000-entry panel (current) |
|---|---|---|---|
| Sampled entries | 100 | 500 | **1,000** |
| Ground-truth mentions | 265 | 1,352 | **2,597** |
| Precision | 97.6% | 97.0% | **98.0%** |
| Recall | 92.8% | 87.4% | **87.8%** |
| F1 | 0.952 | 0.919 | **0.926** |
| Entry-exact | 84.0% (84/100) | 69.6% (348/500) | **71.4% (714/1,000)** |

Recall and entry-exact accuracy both settled well below the original
100-entry panel's numbers once the sample grew large enough to surface
the true long tail, and then held roughly steady between the 500- and
1,000-entry draws (87.4%->87.8% recall, 69.6%->71.4% entry-exact) --
consistent with the panel having already reached a fairly stable,
representative estimate of the extractor's real-world performance by
n=500, with the 1,000-entry panel mainly tightening the confidence in
that estimate rather than revealing a fresh trend. Precision stayed
consistently high (~97-98%) across every panel size, meaning the
extractor rarely invents a compound that isn't there -- its main
failure mode, confirmed at scale, is *missing or mislabeling* long-tail
reagents (buffer-family collapses like "Bis-Tris propane" to
"Bis-Tris", "AMP-PNP" to "AMP", metal-acetate stripping to a bare
"acetate"), not fabricating them. This is the strongest evidence in the
paper for the LLM-assisted fallback described as future work.

- `heldout_sample_1000_raw.json` -- the 1,000 entries (seed=1000),
  unannotated, with the extractor's output attached.
- `heldout_annotations_1000.json` -- the same 1,000 entries with manual
  ground truth and per-entry TP/FP/FN counts (merged from
  `batches_1000/`).
- `heldout_summary_1000.json` -- the aggregate numbers in the table
  above, computed from the file above; this is the result reported in
  the paper's current Table I.
- `batches_1000/` -- the 20 intermediate 50-entry batches
  (`batch_NN_raw.json` / `batch_NN_annotations.json`) the 1,000-entry
  audit was split into.
- `build_sample_1000.py` -- draws the 1,000-entry sample, excluding
  every entry already used in the 100-entry and 500-entry panels.

## 665-entry panel (current Table I result)

Between the 1,000-entry panel above and this one, the dictionary and
parsing logic went through an extended, intensive round of fixes driven
by manually reading entries from the 1,000-entry sample itself --
dozens of missing aliases, several structural regex bugs (an alias
boundary guard that let a short numeric-weight alias swallow a prefix
of a longer one, a reagent-alternation backtracking bug that mis-split
multi-word aliases, European decimal-comma parsing, concentration-unit
ordering, clause-splitting edge cases around parentheses and slash-
joined buffer names, and more), and a handful of pre-existing silent
dict-literal duplicate-key collisions. In effect, the 1,000-entry panel
functioned as a sixth development round, the same way the original
60-entry panel did before the first held-out check existed -- so by the
time this round of fixes was done, it no longer measured held-out
performance honestly.

This panel is the fresh held-out check that follows that tuning round:
every one of the **665 entries still unused** in the 100-protein pool
(2,265 total retained hits, with 1,600 already consumed by the
100+500+1,000-entry panels) -- not a sub-sample, but the *entire*
remainder, since 665 is less than 1,000 and a disjoint 1,000-entry draw
is no longer possible from this pool. This also means the 100-protein
pool is now fully exhausted across all four panels combined: every
retained hit it ever produced has been through manual ground-truth
annotation at least once.

All 665 entries were manually audited with the same single-annotator,
read-the-raw-text-first TP/FP/FN methodology, against the
already-fixed code -- **no code changes made in response to what this
panel found.** The audit was split into 14 batches (thirteen of 50, one
of 15) and annotated independently; `batches_665/` holds every batch's
raw sample and annotations.

| | 1,000-entry panel | 665-entry panel (current) |
|---|---|---|
| Sampled entries | 1,000 | **665** |
| Ground-truth mentions | 2,597 | **1,725** |
| Precision | 98.0% | **99.2%** |
| Recall | 87.8% | **93.6%** |
| F1 | 0.926 | **0.963** |
| Entry-exact | 71.4% (714/1,000) | **84.7% (563/665)** |

Precision, recall, F1 and entry-exact accuracy all improved
substantially over the 1,000-entry panel -- expected, since that
panel's own findings are exactly what drove the intervening fixes, and
this panel draws from the same underlying text distribution. The 13
false positives cluster around a few recurring patterns: a
comma-as-thousands-separator PEG weight ("PEG 20, 000" / "PEG 20.000")
still gets parsed as a much smaller, nonexistent "PEG 20"; a
transition-state-analog complex written as "Mg.ADP.Vanadate" gets
correctly identified but reported as two separate entities rather than
one; and a few single-token substring collisions (ethyl acetate read
as the unrelated salt "acetate", bis-tris propane collapsed to plain
Bis-Tris, a hexanediol isomer mapped to the wrong one). The 111 false
negatives are dominated by two things: depositor typos/misspellings
that a literal-alias dictionary cannot anticipate ("ammonoium
sulphate", "Teis-HCl", "trisodium citrarte", "soduim aceteate", and
similar), and a long, open-ended tail of specific small-molecule
ligands, inhibitors and substrates with no dictionary entry at all
(e.g. vinblastine, NU2058, Ap5A, several opaque pharma compound codes).
A recurring, narrower miss is the no-space "PEG<weight>MME" token
(e.g. "PEG2000MME") appearing without any separator before "MME" --
seen across several batches and not yet covered by the dictionary's
PEG-MME alias patterns. These failure patterns are consistent with the
paper's existing argument for the LLM-assisted fallback as future work:
precision stays very high (the extractor rarely invents a compound
that isn't there), while recall is capped by a long tail of spelling
variants and genuinely novel reagent names no fixed dictionary can
fully enumerate.

- `heldout_sample_665_raw.json` -- all 665 unused entries, unannotated,
  with the current extractor's output attached (re-run against the
  code at sample-build time, not carried over from the original
  pipeline run like earlier panels' raw samples).
- `heldout_annotations_665.json` -- the same 665 entries with manual
  ground truth and per-entry TP/FP/FN counts (merged from
  `batches_665/`).
- `heldout_summary_665.json` -- the aggregate numbers in the table
  above, computed from the file above; this is the result reported in
  the paper's current Table I.
- `batches_665/` -- the 14 intermediate batches (`batch_NN_raw.json` /
  `batch_NN_annotations.json`) the audit was split into.
- `build_sample_665.py` -- draws the 665-entry sample: every entry not
  already used in the 100-, 500-, or 1,000-entry panels.
