# Annotation guidelines: crystallization reagents in PDB records

Thank you for helping. You will read 200 short free-text crystallization
descriptions from the Protein Data Bank (`pdbx_details`) and write down,
for each one, **the list of chemical reagents it mentions**.

You are an independent annotator. Please:

- work only from the text in `annotation_sheet.csv`;
- **do not** run Crystal Explorer or any extraction tool on these entries,
  and do not look at anyone else's annotations;
- you may look up an unfamiliar name (a search engine, PubChem, the PDB
  entry page) to decide whether it is a chemical.

Expect about 1-2 minutes per entry, roughly 4-6 hours in total. It is fine
to split it over several sessions.

## What to write

In the `reagents` column, write each reagent **once**, separated by
semicolons, e.g.

```
ammonium sulfate; Tris; PEG 3350; glycerol
```

- Use a clear standard name. Spelling and capitalization do not matter,
  and abbreviations are fine (`NaCl`, `MPD`, `PEG 3350`), but write what the
  chemical *is*, not a typo (`Imidanzol` -> `imidazole`, `AmSulfate` ->
  `ammonium sulfate`).
- Order does not matter. Concentrations, pH and temperature are **not**
  recorded.
- If an entry contains no reagent at all, write `NONE`.
- If you are unsure about a decision, still make it, and explain in
  `comments`. Those notes are very useful.

## What counts as a reagent

**Count** every named chemical anywhere in the text: reservoir, drop,
protein buffer, additive, cryoprotectant or soaking solution. This
includes:

- precipitants (PEGs, salts, MPD, alcohols, ...), buffers, salts,
  additives, reducing agents (DTT, TCEP, beta-mercaptoethanol), detergents,
  cryoprotectants (glycerol, ethylene glycol), metal salts;
- named ligands, inhibitors, substrates, cofactors and soaking compounds,
  including drug names and compound codes that identify a specific
  molecule (`gefitinib`, `NAD`, `HGN-1310`, `cAMP`);
- a single named reagent product that is one defined substance or
  formulation sold as one reagent (`Tacsimate`, `Jeffamine ED-2001`,
  `pentaerythritol propoxylate 5/4`).

**Do not count:**

- proteins and peptides (the target protein, enzymes, antibodies, protein
  partners, peptides);
- unnamed compounds (`inhibitor`, `compound 4`, `ligand`, `nucleotides`);
- mixtures without listed components: screen names and codes
  (`PEGs II B8`, `JCSG+ A1`, `SM5-E1`), Morpheus mixes (`NPS`,
  `buffer system 2`, `precipitant mix 4`), `CryoProtX`, `tryptone`,
  `dead sea water`. **But** if the text lists the components of a mix,
  count each listed component;
- an acid or base mentioned **only** to adjust the pH (`HCl`, `NaOH`,
  `acetic acid` "to pH 5", "titrated with KOH").

## Counting rules

- **One reagent, counted once per entry.** If the same reagent appears
  twice (e.g. in the reservoir and again in the cryo or soak solution),
  write it once.
- **A buffer written with its counter-ion is one reagent:** `Tris-HCl`,
  `Tris/HCl`, `HEPES-NaOH`, `MES:NaOH`, `sodium HEPES` -> one reagent each
  (`Tris`, `HEPES`, `MES`). `Tris` and `Tris-HCl` in the same entry are the
  same reagent, so write it once.
- **Two salts of the same acid used as one buffer are one reagent:**
  `K2HPO4/KH2PO4`, `Na/K phosphate` -> `potassium phosphate`,
  `sodium/potassium phosphate`.
- **Two different compounds are two reagents, even when written as one
  buffer:** `MES/imidazole` -> `MES; imidazole`; `imidazole-malate` ->
  `imidazole; malate`.
- **Keep chemically different forms apart:** `PEG 5000 MME` is not
  `PEG 5000`; `PEG 400` is not `PEG 4000`; `zinc acetate` is not
  `sodium acetate`; `ATP` is not `AMP-PNP`. Write the specific form named in
  the text. If the text only says `PEG` with no size, write `PEG`.
- **Salts:** write the salt as named (`ammonium sulfate`, `lithium
  chloride`, `CdCl2`). A bare anion written as a reagent (`50 mM acetate`)
  is written as given (`acetate`).

## Examples

| pdbx_details | reagents |
|---|---|
| `0.1 M NACITRAT, PH 4.8, 6% PEG8000, 28% MPD` | `sodium citrate; PEG 8000; MPD` |
| `CADMIUM CHLORIDE, MALATE, HCL, pH 4.9` | `CdCl2; malate` (HCl is the titrant) |
| `0.09 M NPS, Morpheus Buffer System 2 pH 7.5, 42% Morpheus Precipitant Mix 4` | `NONE` |
| `25% PEG 3350, 0.2 M Li2SO4, 0.1 M Bis-Tris/Imidazol pH 6.5, 2 mM DTT. Crystals soaked with compound 7` | `PEG 3350; lithium sulfate; Bis-Tris; imidazole; DTT` |
| `Protein: 8 mg/ml kinase, 2 mM inhibitor, 5% DMSO, 150 mM KCl, 20 mM HEPES-NaOH pH 7.5` | `DMSO; KCl; HEPES` |
| `0.2 M ammonium fluoride, 20% PEG3350 and 4 mM DHP77` | `ammonium fluoride; PEG 3350; DHP77` |

When you are done, send back `annotation_sheet.csv` with the `reagents` and
`comments` columns filled in.
