# 🎉 v1.0.0 — Protein Crystallization Data Extraction (PCDE)

**Protein Crystallization Data Extraction Tool**  
*by Nana Njantang Ruth · ORCID [0000-0002-6003-7521](https://orcid.org/0000-0002-6003-7521)*

---

## Overview

This is the first stable release of the automated pipeline for retrieving, filtering, and analyzing protein crystallization conditions from the RCSB Protein Data Bank (PDB). Starting from an amino acid sequence, the pipeline produces structured, analysis-ready datasets and publication-quality figures.

---

## What's included in this release

### 🗂️ Input & FASTA Export
- Accepts any amino acid, DNA, or RNA sequence interactively
- Automatically detects sequence type (protein / DNA / RNA)
- Saves the input sequence as a FASTA file to `output/{seq_type_name}/`

### 🔍 Step 1 — RCSB Sequence Search (`rcsb_sequence_identity.py`)
- Queries the RCSB Search API (`/rcsbsearch/v2/query`) using sequence similarity scoring
- Maps sequence type to the correct target (`pdb_protein_sequence`, `pdb_dna_sequence`, `pdb_rna_sequence`)
- Filters results to **X-ray crystallography only** via a per-hit lookup to `/rest/v1/core/entry/{pdb_id}`
- Extracts per-hit: `PDB_ID`, `Entity`, `Score`, `Seq_id` (sequence identity), `E-value`
- Saves ranked hits to `{seq_type_name}_rcsb_hits.csv`

### 🧹 Step 2 — Crystallization Data Retrieval & Filtering (`PDB_searchAPI.py`)
- Runs a full PDB search for crystallization metadata
- Filters entries by experimental conditions via `filter_experimental_conditions()`
- All intermediate temporary CSV files are automatically cleaned up after each stage

### 🧪 Compound Annotation (`extract_structures.py`)
- Enriches filtered crystallization data with a `COMPOUND` column using `structures.pkl`
- Output stored in a temporary CSV before merging

### 🔗 Step 3 — Merge RCSB + Crystallization Data
- Merges RCSB sequence-identity results with crystallization condition data on `PDB_ID`
- Final merged CSV contains a defined column order:
  ```
  PDB_ID · Entity · Score · Seq_id · E-value · Resolution · Pubmed_id · Method · pH · Temp · Ligands · Polymer · Assembly · pdbx_pH_range · pdbx_details · Compounds (mM) · PEG_Id · PEG_con
  ```
- Saved as `{seq_type_name}_merged_results.csv`

### 📊 Step 4 — High-Resolution Visualization (`plot.py`)
- `run_plot()` generates **300 dpi** analytical scatter plots from the merged CSV:
  - pH vs. Temperature (K)
  - pH vs. PEG concentration (%)

### 📄 Consolidated PDF Report
- Outputs `Cryst_cocktail_Table.pdf`: a colored summary table compiling PDB IDs, alignment metrics, ligands, and complete chemical cocktails into a single laboratory resource

---

## Pipeline Flow

```
Sequence input
     │
     ├──► FASTA export
     │
     ├──► Step 1: RCSB sequence search + X-ray filter → rcsb_hits.csv
     │
     ├──► Step 2: PDB crystallization search (×6 workers) → temp CSV
     │                └── filter experimental conditions → temp CSV
     │                        └── append COMPOUND column → temp CSV
     │
     ├──► Step 3: Merge RCSB + crystallization data → merged_results.csv
     │
     └──► Step 4: run_plot() → scatter plots (300 dpi)
```

---

## Output Structure

```
output/
└── {seq_type_name}/
    ├── {seq_type_name}_sequence.fasta
    ├── {seq_type_name}_rcsb_hits.csv
    ├── {seq_type_name}_merged_results.csv
    ├── Cryst_cocktail_Table.pdf
    └── plots/
        ├── pH_vs_Temperature.png
        └── pH_vs_PEG_concentration.png
```

---

## Module Overview

| Module | Role |
|--------|------|
| `main.py` | Pipeline orchestration and I/O |
| `rcsb_sequence_identity.py` | RCSB sequence search + X-ray filter |
| `PDB_searchAPI.py` | Crystallization metadata retrieval & filtering |
| `extract_structures.py` | Compound annotation via structures.pkl |
| `plot.py` | High-resolution scatter plot generation |

---

## Dependencies

- Python 3.8+
- `requests`, `pandas`, `openpyxl`, `matplotlib`, `concurrent.futures`, `tempfile`

Install all dependencies with:

```bash
pip install -r requirements.txt
```

---

## 📦 Download Assets

The following pre-built archives are available for this release:

| Asset | Format | Size | Description |
|-------|--------|------|-------------|
| `Source code (zip)` | ZIP | - | Complete repository as ZIP archive |
| `Source code (tar.gz)` | TAR.GZ | - | Complete repository as compressed TAR archive |

**Direct download links:**
- [Download ZIP](https://github.com/RitAreaSciencePark/Protein_Crystallization_Data_Extraction/archive/refs/tags/V1.0.0.zip)
- [Download TAR.GZ](https://github.com/RitAreaSciencePark/Protein_Crystallization_Data_Extraction/archive/refs/tags/V1.0.0.tar.gz)

---

## How to Install

### Option 1: From ZIP Archive
```bash
# Download and extract
unzip V1.0.0.zip
cd Protein_Crystallization_Data_Extraction

# Create virtual environment
python -m venv venv
source venv/bin/activate  # Linux/macOS
venv\Scripts\activate     # Windows

# Install dependencies
pip install -r requirements.txt
```

### Option 2: Clone from Git
```bash
git clone -b V1.0.0 https://github.com/RitAreaSciencePark/Protein_Crystallization_Data_Extraction.git
cd Protein_Crystallization_Data_Extraction
pip install -r requirements.txt
```

---

## Quick Start

### Web Application
```bash
cd protein_crystallization_app
python manage.py migrate
python manage.py runserver
```
Then visit `http://127.0.0.1:8000`

### Command-Line (Single Sequence)
```bash
cd src
python Main.py
Enter sequence: MSPRKTYILKLYVAGNTPNSVRALK...
Enter a descriptive sequence type name: MyProtein
```

### Command-Line (FASTA Batch)
```bash
cd src_fasta_file
python main.py input.fasta
```

---

## How to Cite

If you use this pipeline in your research, please cite:

> Nana Njantang Ruth, Valerio PIOMPONI, and Adrea DALLE VEDOVE (2026). *Protein Crystallization Data Extraction Tool* (v1.0.0). GitHub Repository. https://github.com/RitAreaSciencePark/Protein_Crystallization_Data_Extraction

Or use the `CITATION.cff` file included in this repository.

---

## Known Limitations

- Compound annotation depends on the static `structures.pkl` reference file
- Free-text crystallization fields in PDB are not fully standardized across all entries
- Web interface improvements planned for v2.0.0

---

## License

MIT © 2025 RitAreaSciencePark

This project is released under the MIT License. See the [LICENSE](LICENSE) file for details.

---

## Support & Feedback

For issues, feature requests, or questions:
- [Open an Issue](https://github.com/RitAreaSciencePark/Protein_Crystallization_Data_Extraction/issues)
- [Start a Discussion](https://github.com/RitAreaSciencePark/Protein_Crystallization_Data_Extraction/discussions)

---

**Thank you for using PCDE! 🧬**
