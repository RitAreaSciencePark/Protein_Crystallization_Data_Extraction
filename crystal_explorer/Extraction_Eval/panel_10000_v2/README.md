# Second 10,000-entry evaluation panel (held-out test of the fixed extractor)

The first 10,000-entry panel (`../panel_10000/`) found a set of errors that
were then fixed in `compound_extraction.py` (20 Oct 2026). Because those fixes
were made by looking at that panel, it became a development set. This panel
is a fresh sample built the same way, with no PDB entry and no protein cluster
shared with any earlier panel, to measure the fixed extractor without bias.

## Result (1,000 manually audited entries)

| Metric | Value | 95% CI (bootstrap, 5,000 resamples) | First panel (before fixes) |
|---|---|---|---|
| Precision | **0.993** | 0.989 – 0.996 | 0.992 |
| Recall | **0.942** | 0.934 – 0.951 | 0.927 |
| F1 | **0.967** | 0.962 – 0.972 | 0.958 |
| Entry-level exact match | **84.6 %** (846 / 1,000) | 82.4 – 86.8 % | 81.6 % |

2,898 ground-truth reagent mentions: 2,731 true positives, 20 false
positives, 167 false negatives (`audit_summary.json`).

Automatic coverage over the full 10,000-entry panel (no ground truth needed,
`panel_10000_coverage.json`): 95.1 % of entries yield at least one compound,
28,146 compounds in total (2.81 per entry), 82.1 % of them with a
concentration, 305 distinct compounds.

## How it was built

Same procedure and scripts as `../panel_10000/`, with new seeds:

1. **Targets** (`select_targets.py`, seed 2028): 30 %-sequence-identity
   clusters of X-ray protein entities (60–900 aa) with at least 20 members.
   Excluded: every PDB ID from earlier panels, including all hits retrieved
   by the first 10,000 panel (12,570 IDs), and every cluster used there.
   948 clusters were eligible; 800 were drawn, one representative each.
2. **Pipeline run** (`run_panel.py`): real Crystal Explorer stage 1 + 2
   (identity 0.5, E-value 1.0, 30 hits), 3 shards, 300 s per-target timeout.
   The run was stopped once 683 targets were done, since the pool already
   held more than 10,000 entries.
3. **Panel and audit sample** (`build_samples.py`): pool of 14,145 unique,
   unused entries with `pdbx_details`; 10,000 drawn (seed 2028) =
   `panel_10000_raw.json`; 1,000 of those (seed 2029) = `audit_1000_raw.json`,
   split into `batches/batch_01..20_raw.json`.
4. **Manual audit** (`annotations/batch_NN_annotations.json`): the raw text of
   each entry was read, the true reagents written down and the extractor's
   output scored, with a note for every error. `analyze_audit.py` aggregates.

**Held-out discipline:** the SHA-256 of `compound_extraction.py`
(`c15bbf49…`) was recorded when the samples were built (`panel_meta.json`)
and checked again after the last batch: it is identical, so the extractor
was not changed while this panel was being scored.

## Scoring rules

Identical to the first panel: one TP per distinct reagent per entry; a name
in the right chemical family is correct (generic `PEG`, `acetate` for
uranium acetate, `citrate` for citric acid); a *different* compound is one
FP plus one FN; proteins, peptides, unnamed ligands, named commercial mixes
(Morpheus, CryoProtX, screen codes) and pH titrants are excluded; only
compound identity is scored, not concentration; single annotator.

## Where the errors are

**False negatives (167)**, by kind:

- Ligands, inhibitors and soaks with no dictionary entry (the largest
  group): e.g. HGN-1195/1196/1310, IMP-1002, gefitinib, pelitinib, cAMP and
  Rp/Sp-cAMPS, ddCTP, UDP-Gal, 2-ketoglutarate/2OG, ascorbic acid, ribose,
  maltohexaose, guaiacol, p-vanillin, 1-naphthol.
- Reagents missing from the dictionary: succinic acid, guanidine,
  ethanolamine, K3PO4, YbCl3, sodium oxalate, sodium fumarate, sodium
  phthalate, sodium cholate/glycocholate, deoxycholic acid,
  triethylene glycol, 18-crown-6, gamma-valerolactone, Zwittergent 3-12,
  Jeffamine "JM600", piperazine, PCTP/KMB/AAB buffers.
- Abbreviations and spellings: `CdCl`, `CaCl`, `Na Br`/`NA BR`,
  `Na Thiocyanate`, `Na Acet`, `CH3COONA`, `ACNH4`, `AmSulfate`,
  `Ca(OAc)2`, `MgAcetate`, `NAMES`, `NABICINE`, `NACITRAT`, `HEPS`,
  `IPN`, `PMME5000`, `MEPEG2000`, `PEG 5KMME`, `PEG MME 5K`, `PEG3.4K`,
  `P3350`, `Na-malonate`, `natrium malonate`, `(NH4)H2PO4`, `TEA`.
- Depositor typos and run-together words: `Imidanzol`, `ACTEATE`,
  `AMMONIUMS SULFATE`, `AMMNONIUM`, `Ammonium sulate`, `POTASIUM`,
  `CLORIDE`, `Cacodilate`, `Sodiumformiate`, `Sodiumcitrate`,
  `Lithiumsulfate`, `AMMONIUMSULPHATE`, `dithiothrietol`, `threalose`,
  `2-METHYL-2, 4-PENTANEDIOL`, `HEPES/NaOHpH`.

**False positives (20)** are, again, nearly all a related but different
compound: PEG MME reported as plain PEG when the MME is written in an
unusual position (4: "PEG 5K monomethylether", "MONOMETHYL-PEG 5000", "MME
PEG 2000", "polyethylene glycol 2000 MME"); glycyl-glycyl-glycine read as
glycylglycine (3); a shorter name inside a modified compound
(pyridoxal-5'-phosphate as phosphate buffer ×2, deaza-PLP as PLP,
2-fluoro-cellotriose, myristoyl-CoA, PR-ATP, a Fuc-GlcNAc disaccharide); and
single cases (potassium cacodylate as sodium cacodylate, tert-butanediol as
2,3-butanediol, the long name of bis-tris propane as Tris, "ATP" from "ATP
analog", a misspelled "ploy(ethylene glycol) monomethyl ester 550" as
ethylene glycol, "2OG" as octyl glucoside).

Compared with the first panel, recall rose by 1.5 points and exact match by
3 points while precision stayed at 0.99. The confidence intervals of the two
panels overlap only slightly (recall 0.917–0.937 vs 0.934–0.951; F1
0.952–0.964 vs 0.962–0.972). The two panels are different samples of
proteins, so this is consistent with a real but modest gain, not proof of
one.

## Files

| File | Content |
|---|---|
| `targets.json`, `selection_log.json` | the 800 selected targets and selection parameters |
| `results_shard*.json` | pipeline status per target (the per-target `runs/<target>/` CSVs were deleted on 20 Oct 2026 to save space; `run_panel.py` regenerates them) |
| `panel_10000_raw.json` | the 10,000-entry panel with extractor output |
| `panel_10000_coverage.json` | automatic coverage statistics for the full panel |
| `audit_1000_raw.json`, `batches/` | the audited sample |
| `annotations/` | ground truth and TP/FP/FN per entry, with notes |
| `audit_summary.json` | final metrics, confidence intervals, extractor SHA-256 |
| `panel_meta.json` | seeds, sizes, extractor SHA-256 at build time |
