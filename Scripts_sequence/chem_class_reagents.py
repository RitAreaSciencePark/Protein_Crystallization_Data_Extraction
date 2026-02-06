import pandas as pd
import re

# -------------------------
# Reagent classification
# -------------------------
REAGENT_CLASSES = {
    "buffer": ["hepes", "tris", "mes", "imidazole", "cacodylate", "piperazine"],
    "peg": ["peg", "ethylene glycol", "glycerol"],
    "salt": ["sodium", "potassium", "magnesium", "calcium", "chloride", "acetate", "ammonium", "phosphate", "na", "k", "mg", "ca", "cl"],
    "ligand": ["atp", "bme", "nad", "fad", "cofactor"],
    "solvent": ["water", "ethanol", "methanol", "dioxane"]
}

# -------------------------
# Synonyms / chemical normalization
# -------------------------
SYNONYMS = {
    "na": "sodium",
    "k": "potassium",
    "mg": "magnesium",
    "ca": "calcium",
    "cl": "chloride",
    "bme": "beta-mercaptoethanol",
    "edo": "ethylene glycol",
    "gol": "glycerol",
    "trizma": "tris",
    "naoac": "sodium acetate",
    "naac": "sodium acetate",
    "mgact": "magnesium acetate",
    "mgac": "magnesium acetate",
    "po4": "phosphate"
}

# -------------------------
# Common protein abbreviations
# -------------------------
PROTEIN_ABBREVS = ["kaiB", "kaiC", "bsa", "lysozyme", "myoglobin"]

# -------------------------
# Helper functions
# -------------------------
def normalize_name(name):
    """
    Normalize synonyms to a standard name.
    """
    name = name.lower().strip()
    # Handle acetate patterns Mg(Ac)2 or NaOAc
    name = re.sub(r'mg\s*\(?ac\)?2?', 'mgac', name, flags=re.IGNORECASE)
    name = re.sub(r'na\s*oac', 'naoac', name, flags=re.IGNORECASE)
    return SYNONYMS.get(name, name)

def extract_proteins(text):
    """
    Extract protein names with concentrations.
    Example: "6.5 mg/mL KaiB-KaiC complex, 3 mg/mL BSA"
    """
    proteins = []
    pattern = re.compile(r'(\d*\.?\d+\s*(mg/ml|g/l|µg/ml))\s*([A-Za-z0-9\-]+(?: [A-Za-z0-9\-]+)*)', re.IGNORECASE)
    matches = pattern.findall(text)
    for conc, unit, name in matches:
        proteins.append(f"{name.strip()} ({conc.strip()})")
    
    # Also check for chemical formula proteins
    formula_matches = re.findall(r'\b([A-Z][a-z]?[A-Z]?[a-z]?\d*)\b', text)
    for f in formula_matches:
        if f.lower() in PROTEIN_ABBREVS and f not in proteins:
            proteins.append(f)
    return "; ".join(proteins)

def extract_reagents(text):
    """
    Extract reagents with concentrations and classify into columns.
    Handles chemical formulas and acetate variants.
    """
    reagents_found = {cls: [] for cls in REAGENT_CLASSES}
    
    # Split by commas, 'and', '+'
    parts = re.split(r',| and | \+ ', text, flags=re.IGNORECASE)
    
    for part in parts:
        # Extract concentration if exists
        conc_match = re.search(r'(\d*\.?\d+\s*(mM|M|µM|uM|mg/ml|g/L))', part, re.IGNORECASE)
        conc = conc_match.group(1).strip() if conc_match else ""
        
        # Remove concentration from name
        name = re.sub(r'(\d*\.?\d+\s*(mM|M|µM|uM|mg/ml|g/L))', '', part, flags=re.IGNORECASE).strip()
        name = normalize_name(name)
        
        # Handle chemical formulas (NaCl, Mg(Ac)2)
        formula_match = re.findall(r'([A-Z][a-z]?[A-Z]?[a-z]?\d*)', name)
        if formula_match:
            name = " ".join([normalize_name(f) for f in formula_match])
        
        # Classify reagent
        for cls, keywords in REAGENT_CLASSES.items():
            if any(kw in name.lower() for kw in keywords):
                reagents_found[cls].append(f"{name} ({conc})" if conc else name)
    
    return reagents_found

# -------------------------
# Main processing
# -------------------------
def process_csv(input_csv, output_csv):
    df = pd.read_csv(input_csv)
    
    rows = []
    
    for _, row in df.iterrows():
        pdb_id = row["PDB_ID"]
        details = str(row.get("pdbx_details", ""))
        
        protein_conc = extract_proteins(details)
        reagents_info = extract_reagents(details)
        
        rows.append({
            "PDB_ID": pdb_id,
            "protein_concentration": protein_conc,
            "buffer": "; ".join(reagents_info["buffer"]),
            "peg": "; ".join(reagents_info["peg"]),
            "salt": "; ".join(reagents_info["salt"]),
            "ligand": "; ".join(reagents_info["ligand"]),
            "solvent": "; ".join(reagents_info["solvent"]),
        })
    
    out_df = pd.DataFrame(rows)
    out_df.to_csv(output_csv, index=False)
    print(f"Output CSV saved: {output_csv}")

# -------------------------
# Example run
# -------------------------
if __name__ == "__main__":
    input_csv = "pdb_mmcif_extracted.csv"  # your input CSV
    output_csv = "protein_reagents_parsed.csv"
    process_csv(input_csv, output_csv)
# ----

