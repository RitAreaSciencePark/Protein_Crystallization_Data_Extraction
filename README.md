
<table>
<tr>
<td><img src="crystal_explorer/viewer/static/viewer/crystal-explorer-mark.png" width="90" alt="Project logo"></td>
<td><h1>Protein Crystallization Data Extraction (PCDE)</h1></td>
</tr>
</table>

[![license](https://img.shields.io/badge/license-MIT-green)](LICENSE)
![python](https://img.shields.io/badge/python-3.12-blue)
![django](https://img.shields.io/badge/django-6.0-darkgreen)
![FAIR data](https://img.shields.io/badge/FAIR-data-orange)
![data source](https://img.shields.io/badge/data%20source-RCSB%20PDB-informational)
![plots](https://img.shields.io/badge/plots-Plotly-purple)
![method](https://img.shields.io/badge/method-X--ray%20crystallography-teal)
![Biopython](https://img.shields.io/badge/Biopython-1.87-yellow)
![metadata](https://img.shields.io/badge/metadata-JSON--LD-lightgrey)
![topic](https://img.shields.io/badge/topic-protein%20crystallization-blueviolet)

Protein crystallography remains one of the primary methods for determining three-dimensional protein structures, and identifying the right crystallization conditions is a critical bottleneck. Given a protein, DNA or RNA sequence, this tool finds its X-ray-solved homologs in the RCSB Protein Data Bank (PDB) and turns their free-text crystallization notes into structured reagent and concentration data, with interactive plots, a conditions table and FAIR-compliant metadata. It is available as standalone scripts and as **Crystal Explorer**, a Django web application that runs the whole pipeline from a browser.

---

## Table of Contents

- [Repository Contents](#repository-contents)
- [Project Overview](#project-overview)
- [Key features](#key-features)
- [Workflow Architecture](#-workflow-architecture)
- [How it works](#how-it-works)
- [Project Architecture](#project-architecture)
- [Installation](#installation)
- [Usage](#usage)
  - [Web application](#web-application)
  - [Command-line for standalone pipeline scripts](#command-line-for-standalone-pipeline-scripts)
- [Output files](#output-files)
- [Module reference](#module-reference)
- [Dependencies](#dependencies)
- [Licence](#licence)
- [Citation](#citation)

---

## Repository Contents

### `Input/`

The scripts import lists and dictionaries from files in **json** format using the json Python module. The description of the **json** files in this directory are described [here](https://github.com/maxdudek/crystallizationDatabase/wiki/Explanation-of-Files)

### `src/`

This directory contains all the python codes for the general pipeline handling a single sequence as input.

### `src_fasta_file/`

This directory contains the codes that handle multi sequences fasta file.

### `Structures/`

This directory contains **structure.pkl** file which is the pdb database with the list of compound extracted from the **pdbx_details** at the level of the **exptl_crystal_grow** section of the mmCIF file. 

### `crystal_explorer/`

This directory contains **Crystal Explorer**, the Django web application: a form-driven front end that runs a four-stage pipeline (RCSB search → compound extraction → grouping and plotting → FAIR metadata) in the background and shows the results as interactive plots and a conditions table. See [`crystal_explorer/README.md`](crystal_explorer/README.md) for setup and project-layout details specific to the app.

### `crystal_explorer/Extraction_Eval/`

Evaluation of the web app's compound extractor against manually audited PDB entries.


--- 
## Project Overview

Determining the right conditions to crystallize a protein is one of the most time-consuming challenge in structural biology. The primary goal is to eliminate the manual bottleneck of screening thousands of chemical cocktails by providing a **PDB-assisted shortcut**. This project addresses that challenge by mining the PDB for crystallization conditions used to solve structures of proteins similar to a given query sequence, allowing researchers to skip initial trial-and-error phases and move directly to optimization.

---
## Key features

### Automated sequence-based PDB mining
The pipeline accepts any protein, DNA, or RNA sequence and automatically queries the RCSB PDB using the MMseqs2 sequence similarity engine. It filters results to X-ray diffraction structures only and ranks hits by sequence similarity score, eliminating manual database browsing entirely.

### Two-tier metadata extraction
Every crystallization parameter is extracted using a two-tier strategy: the structured mmCIF field is read first, and if absent or empty, regular expression mining is applied to the free-text `pdbx_details` field. This ensures that no experimentally reported information is lost due to incomplete structured annotation.

### Intelligent free-text parsing
The free-text `pdbx_details` field is turned into compound–concentration pairs. The standalone scripts (`src/`) use the crystallization database of Lynch, Dudek et al. (`structures.pkl`); the web app uses its own curated dictionary (`compound_extraction.py`) of **409 reagents and 1,629 aliases** (abbreviations, formulas, common typos), each classified as Buffer, Salt, Precipitant or Additive.

### Parallel data collection
Metadata extraction uses a configurable multi-threaded `ThreadPoolExecutor`, processing multiple PDB entries simultaneously. A live progress bar (tqdm) updates in real time as each entry completes, giving immediate visibility into pipeline progress.

### FAIR-compliant data provenance
Every row in the output dataset is traceable to its source through the PDB accession code, PubMed ID, and experimental method. All data originates from the publicly accessible RCSB REST and GraphQL APIs using open standards (mmCIF, JSON, FASTA), and outputs are saved in non-proprietary formats (CSV, PDF, PNG).

### Interactive web application
**Crystal Explorer**, a Django web application, wraps the full pipeline in a browser interface: background searches, interactive Plotly plots and reagent-frequency charts, a conditions table linked to RCSB, UniProt, PubChem and a Mol\* 3D view, a per-session search history, and downloadable outputs (PDF, PNG, CSV, FAIR JSON-LD metadata).

### High - Resolution Vizualization
Two scatter plots are generated automatically: pH versus temperature and pH versus PEG concentration with data points coloured by sequence similarity score using the viridis colormap and shaped by crystallization method. A coloured PDF summary table compiles successful crystallization conditions with full metadata.

## 🧬 Workflow Architecture

The pipeline is organized into four main stages, as illustrated in the project’s graphical abstract:

### 🔹 Step 1: Sequence Input & PDB Search

- **Input**  
  Accepts a raw protein sequence (plain string or FASTA format).

- **Categorization**  
  The `detect_seq_type` function determines whether the input is **protein, DNA, or RNA**.

- **Search**  
  Sends an automated query to the **RCSB PDB Search API** using identity and e-value thresholds to identify homologous structures solved via **X-ray diffraction**.

---

### 🔹 Step 2: Data Extraction & Filtering

- **Parallel Processing**  
  Utilizes six parallel workers to efficiently handle large-scale data retrieval.

- **Caching**  
  Implements local caching for **mmCIF files** and **PubMed IDs** to improve performance.

- **Extraction Logic**  
  The `extract_mmcif_info` module aggregates key metadata, including:
  - Resolution  
  - Polymer type (`uni_pol` vs. `complex`)  
  - Assembly details (monomer, dimer, etc.) etc.

- **Filtering**  
  The `filter_experimental_conditions` function removes entries missing critical experimental parameters which are **pH, temperature, method and pdbx_details**.

---

### 🔹 Step 3: Chemical Enrichment and Normalization

- **Compound Mapping**  
Maps PDB IDs to a specialized dictionary (structures.pkl) and extract the list of compounds from the dictionary to enrich the csv file. From the list of compounds, PEG_ID and the PEG_concentration were extrated for futher analysis.

-**Publication Grouping**
 PDB ID having the same Conditions (PubMed ID, pH, temperatures and compounds) are grouped in one line. 

---

### 🔹 Step 4: Plotting & Reporting

The `run_plot()` stage handles data visualization and output generation:

- **Visual Encoding**  
  - Crystallization methods → unique markers  
  - Alignment scores → `viridis` colormap  

- **Generated Plot**
  - pH vs. PEG concentration (%), including a *"No PEG"* and *"No pH"* categories. 

- **Final Report**  
  Produces a color-coded PDF table: {name}_Cryst_cocktail_Table.pdf groups unique conditions by publication and experimental parameters

---

## How it works

```
Query sequence (FASTA or plain string)
        │
        ▼
RCSB MMseqs2 sequence search
(identity cutoff 30 %, X-ray only)
        │
        ▼
Retrieve crystallization metadata
(REST API + mmCIF parsing via gemmi)
        │
        ▼
Filter · merge · standardize
(pandas · compound parsing · PEG extraction)
        │
        ▼
Vizualisation
(pH vs PEG scatter · PDF cocktail table)
        │
        ▼
Output files
(CSV · FASTA · PNG plots · PDF table)
```

The web application (`crystal_explorer/`) wraps the entire pipeline in a Django interface. The search view hands the pipeline to a background thread and redirects to the results page immediately; that page shows a "Search running..." panel that auto-refreshes until the run finishes (or shows the error if it failed).

---


## Project Architecture

```text
┌──────────────────────────────────────────────────────────┐
│       Protein Crystallization Data Extraction Tool       │
└────────────────────────────┬─────────────────────────────┘
                             │
         ┌───────────────────┼───────────────────┐
         │                   │                   │
   ┌─────▼─────┐       ┌─────▼─────┐       ┌─────▼─────┐
   │   src/    │       │ src_fasta │       │  Django   │
   │   Core    │       │   File    │       │  Web App  │
   │  Scripts  │       │  Handler  │       │           │
   └─────┬─────┘       └───────────┘       └─────┬─────┘
         │                                       │
         ├── PDB Search API                      ├── Web Interface
         ├── Data Extraction                     ├── Pipeline Views
         ├── Metadata Generation                 ├── URL Routing
         ├── Database Operations                 └── Forms & Utils
         ├── Sequence Analysis                         │
         └── Visualization                       ┌─────▼─────┐
                                                 │  SQLite   │
                                                 │ Database  │
                                                 └───────────┘
```

---

## Installation

### Prerequisites

- Python 3.12 or later
- pip

### Steps

```bash
# 1. Clone the repository
git clone https://github.com/RitAreaSciencePark/Protein_Crystallization_Data_Extraction.git
cd Protein_Crystallization_Data_Extraction

# 2. Create and activate a virtual environment
python -m venv venv
source venv/bin/activate          # Linux / macOS
venv\Scripts\activate             # Windows

# 3. Install dependencies
pip install -r requirements.txt                    # standalone scripts (src/, src_fasta_file/)
pip install -r crystal_explorer/requirements.txt   # web application

# 4. Apply database migrations (web app only)
cd crystal_explorer
python manage.py migrate

# 5. Start the development server
python manage.py runserver
```

The web application will be available at `http://127.0.0.1:8000`.

---

## Usage

### Web application

1. Open `http://127.0.0.1:8000` and go to **Explorer**.
2. Enter a **protein name** (used to name the output folder/files) and paste your **sequence** (raw or FASTA).
3. Choose the **sequence type** (protein / DNA / RNA) and, optionally, adjust the identity, E-value, and max-hits thresholds.
4. Click **Run search**.
5. The results page shows a progress panel until the background search finishes.
6. Browse the interactive plots, reagent-frequency charts and conditions table, or download the static PDF/PNG/CSV/metadata exports. Past runs are listed under **History** (per browser session) and are reused automatically if you resubmit the same sequence and settings.

---

### Command-line for standalone pipeline scripts

#### Single sequence

```bash
cd src
python Main.py
Enter sequence: MSPRKTYILKLYVAGNTPNSVRALK...
Enter a descriptive sequence type name: KaiB
```
---

#### FASTA file input 

```bash
cd src_fasta_file
python main.py "NAME_OF_THE_FASTA_FILE".fasta
```
Produces one CSV file per sequence in the FASTA file. Each CSV contains the PDB ID, similarity score, and experimental crystallization data for all homologous structures found.

---
#### FASTA file input — combined CSV (no duplicates)

```bash
cd src_fasta_file
python main2.py "NAME_OF_THE_FASTA_FILE".fasta
```
Produces a single combined CSV file for all sequences in the FASTA file. PDB entries that appear as hits for multiple sequences are included only once, making this output suitable for batch analysis across a protein family.

---

## Output files

| File | Format | Description |
|---|---|---|
| `{name}_sequence.fasta` | FASTA | Input sequence saved in standard format |
| `{name}_rcsb_hits.csv` | CSV | MMseqs2 search results: PDB ID, entity, score, identity (%), E-value|
| `{name}_merged_results.csv` | CSV | Full merged dataset: crystallization conditions + sequence similarity scores |
| `{name}_PEG.png` | PNG | Scatter plot: pH vs PEG concentration (%), coloured by similarity score |
| `{name}_Cryst_cocktail_Table.pdf` | PDF | Coloured summary table of all unique crystallization conditions |

### Merged CSV columns

| Column | Description |
|---|---|
| `PDB_ID`| The unique 4-character PDB accession code |
| `Entity` | Entity number within the PDB entry |
| `Score` | RCSB MMseqs2 similarity score (0–1) |
| `Seq_id` | Sequence identity (%) |
| `E-value` | Statistical significance of the alignment |
| `Resolution` | Diffraction resolution (Å) |
| `Pubmed_id` | PubMed identifier of the associated publication |
| `Method` | Crystallization method (e.g. sitting drop vapour diffusion) |
| `pH` | Crystallization pH |
| `Temp` | Crystallization temperature (K) |
| `Ligands` | Co-crystallized ligands, cofactors, or inhibitors |
| `Polymer` | uni_pol (single polymer) or complex (multiple chains) |
| `Assembly` | Oligomeric state from _pdbx_struct_assembly (monomer, dimer, etc.) |
| `pdbx_pH_range` | pH range when a single value is not reported |
| `pdbx_details` | Free-text crystallization details (REMARK 280 equivalent) |
| `Compounds(con_unit=mM)` | Parsed compound list with concentrations in mM |
| `PEG_Id` | PEG molecular weight identifier (e.g. 4000) |
| `PEG_con` | PEG concentration (%) |

---

## Module reference

### `rcsb_sequence_identity.py`

Performs a sequence similarity search against the RCSB PDB using the MMseqs2-based Search API. Filters results to X-ray diffraction structures only. Returns PDB ID, entity, similarity score, sequence identity, and E-value for each hit.

**Key function:** `run_and_save(sequence, output_csv_1)`

---

### `PDB_searchAPI.py`

This module is the primary data collection engine of the pipeline. It queries the RCSB PDB for structures similar to a query sequence, downloads their experimental metadata, and saves the results to a CSV file.

**Sequence search.** `search_pdb_by_sequence()` submits a combined query to the RCSB Search API: a sequence similarity search (MMseqs2, identity cutoff 50 %, E-value ≤ 1×10⁻⁵) intersected with an X-ray diffraction filter, so only crystallographically determined structures are returned. Results are ranked by score and filtered against a configurable score cutoff.

**Parallel metadata extraction.** For each hit, `extract_mmcif_info()` is called concurrently across a thread pool. It downloads and caches the mmCIF file for the entry, then extracts: resolution, PubMed ID, polymer type (single chain vs complex), biological assembly, crystallization method, pH, temperature, free-text crystallization details (`pdbx_details`), pH range, and ligands.

**Two-tier field extraction.** Four dedicated functions — `get_ph_from_mmcif_or_details()`, `get_method_from_mmcif_or_details()`, `get_temperature_from_mmcif_or_details()`, and `get_pdbx_ph_range_from_mmcif_or_details()` — each attempt to read the value from the structured mmCIF field first, and fall back to regular expression mining of the free-text `pdbx_details` string if the structured field is absent or empty.

**Caching.** Downloaded mmCIF files and PubMed ID lookups are cached to a local `.pdb_cache/` directory so repeated runs do not re-download data already on disk.

**Filtering.** `filter_experimental_conditions()` post-processes the output CSV to retain only rows that have at least one of pH, temperature, method, or pdbx_details populated — discarding entries with no usable crystallization information.

**Key functions:**
- `search_pdb_by_sequence(sequence, output_csv, max_workers)`
- `filter_experimental_conditions(input_csv, output_csv)`

---

### `extract_structures.py`

 Lynch et al. developed a Python-based tool in 2020 for creating a searchable and updatable database of crystallization conditions extracted from the free-text crystallization details available in PDB entries. The database includes a parsing method that converts crystallization descriptions into structured lists of chemical compounds and their concentrations, along with a dictionary for standardizing compound names. The tool is freely available through the [GitHub repository](https://github.com/maxdudek/crystallizationDatabase). An updated version of the database was generated for this study by executing the standalone `pdb_crystal_database.py` script located in the `src` subdirectory. The workflow was designed to parse the inconsistently formatted crystallization details associated with PDB structures solved using the X-ray diffraction (XRD) method into a standardized list of compounds suitable for computational analysis. In addition, the pipeline was executed using a newer version of Natural Language Toolkit (NLTK) with the `tokenizers/punkt_tab` tokenizer, enabling the successful extraction of crystallization data from 195,183 structures.
  The database stored in the `structures.pkl` file of the `structures` sub-folder. The database can be updated as the PDB keep on growing in number of structures. 

  The extract_structure.py script matches entries using PDB identifiers and appends standardized compound and concentration information from the structures.pkl database to the filtered CSV dataset. In addition, the workflow automatically extracts PEG-related information, including PEG molecular weight identifiers and concentrations, through regular-expression-based parsing. The final processed dataset is exported as a CSV file containing structured crystallization conditions suitable for downstream computational and statistical analysis.

Appends compound and PEG data from a pre-built compound library (`structures.pkl`) to the filtered crystallization CSV. Uses a fake module injection pattern to safely unpickle `Structure` objects without requiring the original package.

**Key functions:**
- `format_compounds(compound_list)` — converts flat alternating list to `"Compound (concentration)"` string
- `extract_peg_info(compound_str)` — extracts PEG molecular weight and numeric concentration
- `append_compound_to_filtered_csv(structures_file, filtered_csv, output_csv)`

---
### `plot.py`

Generates all visualisation outputs from the merged CSV. Handles missing pH and temperature values using sentinel coordinates, assigns method-specific marker shapes, and colours data points by similarity score using the viridis colormap.

**Key function:** `run_plot(output_csv_file, protein_name)`

**Outputs:**
- `{name}_TEMP.png` — pH vs temperature scatter plot
- `{name}_PEG.png` — pH vs PEG concentration scatter plot
- `{name}_Cryst_cocktail_Table.pdf` — coloured condition summary table

---

### `src_fasta_file/cleaning_and_read_fasta_file.py`

Reads a raw FASTA file containing one or more sequences, removes gap characters (`-`), and returns a dictionary mapping sequence IDs to cleaned sequences.

---

### `src_fasta_file/extract_data_fasta.py`

Performs a PDB sequence search for each sequence in a cleaned FASTA dictionary. Parses the resulting mmCIF files using `gemmi` to extract experimental crystallization details. The identity cutoff is set to 30 % to maximise the number of retrieved structures.

**Key function:** `extract_crystallization(pdb_id)` — pulls all crystallization fields from the mmCIF block via `doc.sole_block()`

---
### For the web application (`crystal_explorer/`)

The web app re-runs the same pipeline logic as its own importable `pipeline/` package (`pdb_sequence_search.py`, `compound_extraction.py`, `plot.py`, plus `metadata_generator.py` for FAIR-compliant output metadata), rather than importing directly from `src/`, so it can evolve independently of the standalone CLI scripts.

#### `pipeline/compound_extraction.py`

Extracts reagents and their concentrations from `pdbx_details` by splitting the text into clauses and matching each against the curated reagent dictionary.

**Key function:** `process_csv(input_csv, output_csv)`

#### `viewer/forms.py`

Defines `SequenceSearchForm`: protein name, sequence, sequence type, identity/E-value/max-hits thresholds, and an optional self-reported lab/institution name.

#### `viewer/views.py`

Handles HTTP requests. `index` spawns the pipeline on a background thread and redirects to the results page immediately (reusing a prior run instead of re-running it if the same sequence/settings were already searched); `results` shows a running/failed state or loads a completed run's `Grouped_conditions.csv` and renders its plots/table; `history` lists past runs; `download` serves the static exports (PDF, PNG, FASTA, compounds CSV, metadata JSON).

#### `viewer/data.py`

Turns `Grouped_conditions.csv` into the interactive Plotly figures and the HTML table rows shown on the results page.

See [`crystal_explorer/README.md`](crystal_explorer/README.md) for the full project layout.

---

## Dependencies

| Package | Purpose |
|---|---|
| `django` | Web application framework |
| `requests` | HTTP calls to RCSB APIs |
| `pandas` | Data manipulation and CSV handling |
| `gemmi` | mmCIF file parsing |
| `matplotlib` | Scatter plots and PDF table generation |
| `numpy` | Numerical operations |
| `tabulate` | Terminal table formatting |
| `pdf2image` | to open the pdf file in the web application |
| `plotly` | Interactive plots in the web application |
| `biopython` | Sequence handling |

Install all dependencies with:

```bash
pip install -r requirements.txt                    # standalone scripts
pip install -r crystal_explorer/requirements.txt   # web application
```
---
## Licence

This project is released under the MIT licence included in the [`LICENSE`](LICENSE) file: `Copyright (c) 2025 RitAreaSciencePark`

---

## Citation

If you use this tool in your research, please cite the associated thesis and the Protein Data Bank:
> R. N. NANA<sup>1</sup>, Valerio PIOMPONI<sup>1</sup> and Adrea DALLE VEDOVE<sup>2</sup> (2026) "Data Extraction Tool for Protein
Crystallization Conditions".

> M. L. Lynch, M. F. Dudek, and S. E. Bowman (2020) A searchable database of crystallization cocktails in the pdb: analyzing the chemical condition
space. *Patterns 1* no. **4**,  







