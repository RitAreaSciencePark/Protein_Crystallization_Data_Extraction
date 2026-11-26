# Retrospective check: do homolog conditions recover known ones?

Paper Section VI-D (`\label{sec:retro}`), Table `tab:retro`, Fig. `fig:retro`.
Written in response to the reviewer comment that the paper did not say
whether suggested cocktails reproduce known conditions.

Uses the held-out 10,000-entry panel (`../panel_10000_v2/`): 674 targets
from distinct 30%-identity clusters, so entries of one target are homologs
and entries of different targets are unrelated proteins.

## Reproduce

```bash
cd crystal_explorer
.venv/bin/python3 Extraction_Eval/retrospective/fetch_metadata.py        # -> entry_metadata.json (RCSB GraphQL, ~2 min)
.venv/bin/python3 Extraction_Eval/retrospective/analyze_retrospective.py # -> retrospective_summary.json, pairs.csv.gz, loo_*.csv
.venv/bin/python3 Extraction_Eval/retrospective/fetch_context.py         # -> entry_context.json (assembly, chains, ligands)
.venv/bin/python3 Extraction_Eval/retrospective/context_analysis.py      # -> context_summary.json, pairs_context.csv.gz
.venv/bin/python3 Extraction_Eval/retrospective/plot_retrospective.py    # -> paper/figures/retrospective.pdf/.png
.venv/bin/python3 Extraction_Eval/retrospective/panel_coverage.py        # -> panel_coverage_stats.json (paper VI-A "Panel composition")
.venv/bin/python3 Extraction_Eval/retrospective/robustness.py            # -> robustness_summary.json, timesplit_*.csv (paper VI-D "Robustness")
```

`panel_coverage.py` also needs `entry_metadata_extra.json` (authors,
crystallization method and source organisms per entry), which was fetched
from the same RCSB GraphQL endpoint with a one-off query (fields
`audit_author`, `exptl_crystal_grow.method`,
`polymer_entities.rcsb_entity_source_organism`); it is kept alongside
`entry_metadata.json`.

Seed 2030; 1,000 cluster-bootstrap resamples over targets. Note that
`entry_metadata.json` is a snapshot of the RCSB annotations at fetch time
(Nov 2026); refetching later may change UniProt mappings slightly.

## Method (short)

- **Entry's protein** = the polypeptide chain with the highest identity
  to the target sequence.
- **Pairwise identity** = global alignment (BLOSUM62, gap -10/-0.5, free
  end gaps), identities / length of the shorter chain.
- **Condition similarity** from the extractor's output: reagent-set Jaccard
  (after merging spelling variants such as Tris-HCl/Tris); shared
  precipitant (PEG of the same MW class <1k / 1-5k / >=6k, an organic
  precipitant such as MPD, or any salt at >= 0.8 M); |delta pH| (mmCIF
  `exptl_crystal_grow.pH`).
- **Baseline** = 180,522 random pairs of entries from different targets.
- **Leave-one-out**: hide an entry and take its top-5 groupmates by identity
  as suggestions. Scenarios: all homologs / other publications only /
  different UniProt only. Baselines: 5 random unrelated entries, and the
  panel's most frequent reagents (size-matched).

## Headline numbers

| | Jaccard | shared precipitant | mean abs dpH |
|---|---|---|---|
| unrelated pairs | 0.057 | 0.294 | 1.44 |
| homolog pairs, 50-70 % id | 0.071 | 0.345 | 1.40 |
| homolog pairs, 90-98 % id | 0.147 | 0.432 | 1.15 |
| homolog pairs, >= 98 % id | 0.474 | 0.695 | 0.64 |
| same UniProt | 0.450 | 0.678 | 0.67 |
| same UniProt, different paper | 0.393 | 0.645 | 0.75 |
| different UniProt | 0.107 | 0.381 | 1.31 |

Spearman identity vs Jaccard: 0.49 [0.46, 0.52] all pairs; 0.16 [0.10, 0.22]
for pairs with different UniProt accessions.

Leave-one-out, precipitant recovered in top-5 (homologs / generic / random):
all homologs 0.864 / 0.658 / 0.646; other publications 0.798 / 0.660 / 0.646;
different UniProt 0.564 / 0.669 / 0.586.

### Robustness (`robustness.py`)

- **Per-target weighting** (mean of per-target means): same conclusions.
  Precipitant recovered 0.828 / 0.754 / 0.543 for the three scenarios,
  against a generic list at 0.641-0.648.
- **Time split** (queries deposited >= 2020-01-01; only earlier entries can
  be suggested; generic list built from earlier entries only):
  - All homologs (2,386 queries, 409 targets): precipitant 0.734 vs generic
    0.686 (paired +4.8 pts [0.7, 8.8]); reagent recall 0.507 vs 0.413
    (+9.4 [5.8, 13.5]); median abs dpH 0.65 vs 1.0.
  - The gain is entirely in the 77% of queries whose closest earlier
    homolog has the same UniProt accession.
  - Different UniProt: generic wins (precipitant 0.573 vs 0.733,
    -16.1 pts [-21.6, -10.5]).

Pairs at >= 98 % identity share a UniProt accession 91-96 % of the time,
69 % at 95-98 % and 49 % at 90-95 %.

**Reading:** strong recovery when the same protein was crystallized before,
including by other groups; no advantage over a frequency prior for a
different protein. This is consistent with Abrahams & Newman (2019).
Recovering deposited conditions does not measure success rate, because the
PDB holds only successes.

### Structural context (`context_analysis.py`)

Written in response to the comment that sequence identity does not imply the
same crystallized object (alone vs. complex, oligomer, construct, ligand).
Context per entry: oligomeric state of assembly 1 (RCSB global symmetry),
the other polymer entities, ligands (non-polymer components minus the
`ADDITIVES` list of ions, buffers, cryoprotectants and glycans); per pair:
alignment coverage = aligned columns / longer chain.

- **Same composition** (both single-entity, same oligomeric state,
  coverage >= 90 %; 51 % of pairs): Jaccard by bin 0.100 / 0.062 / 0.099 /
  0.195 / 0.569; Spearman 0.52 [0.48, 0.56]. The identity trend is not a
  context artifact.
- **Same UniProt, same vs. different context** (partners, oligomeric state,
  coverage): Jaccard 0.565 vs 0.232 (all papers), 0.503 vs 0.218 (different
  papers). One factor differing (different papers): coverage 0.204,
  partners 0.251, oligomeric state 0.346, apo vs ligand 0.369.
- Same context with *different* ligands (0.547) is above same ligands
  (0.456): ligand-soaking series reuse one crystal form.
