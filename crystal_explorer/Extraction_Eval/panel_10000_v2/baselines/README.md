# Baseline: Lynch et al. (2020) parser on the same entries

Reviewer request: compare against the Lynch et al. database/parser on the
same entries.

`run_lynch2020.py` runs the **published, unmodified** parser of Lynch,
Dudek & Bowman (2020, *Patterns* 1:100024) on exactly the same `pdbx_details`
texts as our 1,000-entry held-out audit:

- code: github.com/Hauptman-Woodward/crystallizationDatabase, pinned commit
  `42be6d1171ba1c2487e8cc4d22d86c7efd3bf84c` (8 Jan 2021);
- database: Zenodo doi 10.5281/zenodo.3931013 (December 2019 snapshot).

We use the code rather than the released database because many audited
entries were deposited after December 2019, and because running both tools on
the same text removes any difference in input.

**Not** used: the modified copy in this repository's `src/` and `Input/`
(adapted for this project), nor `Structures/structures.pkl` built with it.
The baseline must be the authors' own version.

## Output

| System | Content |
|---|---|
| `lynch2020_all/` | every compound name their parser returns, recognized or not |
| `lynch2020_recognized/` | only names in their compound dictionary, i.e. what their curated database's compound matrix contains |

Each row also carries their `sensible` flag (every extracted name recognized
and at least one compound). Only sensible entries enter their curated
database.

On the 1,000 audited entries: no parse errors, 75.5 % sensible (their paper:
about 74 %, 99k of 134k), 4.2 % empty, and 3.00 names per entry (2.66
recognized). Details are in `run_meta.json`.

## Reproduce

```bash
git clone https://github.com/Hauptman-Woodward/crystallizationDatabase lynch
git -C lynch checkout 42be6d1171ba1c2487e8cc4d22d86c7efd3bf84c
python3 -m venv lynch_venv && lynch_venv/bin/pip install nltk
export NLTK_DATA=$PWD/nltk_data
lynch_venv/bin/python -c "import nltk; nltk.download('punkt', download_dir='$NLTK_DATA'); nltk.download('punkt_tab', download_dir='$NLTK_DATA')"
lynch_venv/bin/python run_lynch2020.py --lynch-dir lynch --entries all
```

The script refuses to run if the clone isn't at the pinned commit or has
local changes. Their dictionary-maintenance functions rewrite their JSON
files; that file writing is disabled, and the in-memory updates are applied
exactly as in their `standardizeAllNames`.

## Result against annotator A's ground truth (all 1,000 audit entries)

`compare_lynch_vs_A.py` scores the parser against A's audit ground truth,
reconstructed as reagent lists in `../annotator_A_gt/`, next to our
extractor's audited counts. Name equivalences are in
`lynch_name_matches.csv`; they apply the same rules A applied to our
extractor (same reagent or right family counts, a different compound does
not). The reconstruction and the equivalences were written by Claude from A's
notes, so **annotator A must verify both before these numbers are used**.

| (2,899 GT reagents) | Precision | Recall | F1 | Entry exact |
|---|---|---|---|---|
| Our extractor (audit) | 0.993 [0.989-0.996] | 0.942 [0.934-0.950] | 0.967 | 0.846 |
| Lynch parser, every returned name | 0.941 [0.926-0.956] | 0.955 [0.943-0.967] | 0.948 | 0.878 |
| Lynch parser, recognized names only | 0.995 [0.991-0.998] | 0.892 [0.877-0.907] | 0.941 | 0.793 |
| Lynch database view (sensible entries) | 0.998 [0.994-1.000] | 0.708 [0.675-0.741] | 0.828 | 0.773 |

How to read it:

- **Splitting the text vs. identifying the chemical.** Counting every
  string it returns, the parser has higher recall than ours, but 182 of its
  hits are raw text spans its own dictionary doesn't recognize. These are
  often the typos and abbreviations ours misses ("ammoniums sulfate",
  "heps", "cdcl"), counted as correct only because the span names the
  reagent. Its 173 false positives are prose and sample codes (`puck`,
  `mg/m`, `set-up`, crystal IDs).
- **Standardized output.** Once names are standardized, which is what makes
  the output usable, precision matches ours (0.995 vs. 0.993) and recall is
  lower (0.892 vs. 0.942).
- **The database itself** keeps only "sensible" entries (755 of 1,000), so
  as a lookup it covers 0.708 of the reagents. On those 755 entries both
  tools are essentially equal (F1 0.987 vs. 0.988).
- **Scope.** On the 61 entries whose text contains "cryo", "soak" or
  "reservoir solution", its recall falls to 0.71 (ours 0.90), as expected
  from its design. On the other 939 entries, its every-name recall is 0.981.
- **Complementary errors.** Of the 2,899 reagents, both tools find 2,618,
  only ours finds 114, and only the parser finds 151.

Per-entry detail is in `lynch_vs_A_entries.csv`; all numbers are in
`lynch_vs_A_summary.json`.

## Scoring on the blinded subset

`../second_annotator/score_agreement.py` picks up both systems
automatically (`baseline:lynch2020_all`, `baseline:lynch2020_recognized`).
It scores them against annotator B and against the adjudicated consensus on
the 200 blinded entries, with the same matching and adjudication as our
extractor and the LLMs.

During adjudication, use the `note` column to tag each of their misses that
comes from a deliberate scope choice rather than a parsing failure. Their
parser **by design** cuts the text at "cryo" and "soaked", and keeps only the
"reservoir solution" part when that phrase is present. So it doesn't report
cryoprotectants, soaked ligands or protein-buffer components, which our rules
count. Suggested tags: `scope:cryo`, `scope:soak`, `scope:protein-buffer`.
The paper should report recall with and without these.

Also report their coverage: the fraction of entries their curated database
would contain at all (`sensible`), next to precision and recall.
