import pandas as pd
import re
import os

# ---------------- Robust Reagent Extractor ----------------
def extract_reagents(pdbx_details):
    """
    Extract Protein, Salt, Buffer, PEG from free-text pdbx_details.
    Returns a dictionary with lists of strings: 'conc reagent'
    """
    if not pdbx_details or pdbx_details in {".", "?", ""}:
        return {"protein": [], "salt": [], "buffer": [], "peg": []}

    text = pdbx_details.replace("\n", " ").strip().lower()

    # Keyword lists
    buffers = ["tris", "hepes", "mes", "mops", "bis-tris", "cacodylate", "acetate", "phosphate", "citrate"]
    salts = ["nacl", "kcl", "mgcl2", "cacl2", "zncl2", "ammonium sulfate", "lithium sulfate", "sodium sulfate", "potassium chloride"]
    pegs = ["peg", "polyethylene glycol"]

    # Patterns for concentrations (number + unit)
    conc_pattern = r"(\d+(?:\.\d+)?\s*(?:m|mm|mg/ml|%|g/l))"

    # Initialize lists
    proteins_list = []
    salts_list = []
    buffers_list = []
    pegs_list = []

    # Split by common delimiters
    parts = re.split(r"[;,\.]", pdbx_details)

    for part in parts:
        part_clean = part.strip()
        if not part_clean:
            continue

        part_lower = part_clean.lower()

        # --- PEG ---
        for peg_kw in pegs:
            if peg_kw in part_lower:
                conc_match = re.search(r"(\d+(?:\.\d+)?\s*%)", part_clean)
                conc = conc_match.group(1) if conc_match else ""
                pegs_list.append(f"{conc} {part_clean}".strip())
                break  # avoid double-counting

        # --- Buffers ---
        for buf in buffers:
            if buf in part_lower:
                conc_match = re.search(conc_pattern, part_clean, re.I)
                conc = conc_match.group(1) if conc_match else ""
                buffers_list.append(f"{conc} {part_clean}".strip())
                break

        # --- Salts ---
        for salt in salts:
            if salt in part_lower:
                conc_match = re.search(conc_pattern, part_clean, re.I)
                conc = conc_match.group(1) if conc_match else ""
                salts_list.append(f"{conc} {part_clean}".strip())
                break

        # --- Protein (mg/mL or identifier) ---
        if "mg/ml" in part_lower or re.search(r"\b[A-Z]{1,3}[a-zA-Z0-9\.\-]+\b", part_clean):
            proteins_list.append(part_clean.strip())

    return {
        "protein": proteins_list,
        "salt": salts_list,
        "buffer": buffers_list,
        "peg": pegs_list
    }


# ---------------- Read CSV from terminal ----------------
input_csv = input("Enter the input CSV filename (with .csv extension): ").strip()
if not os.path.isfile(input_csv):
    raise FileNotFoundError(f"File '{input_csv}' not found.")

df = pd.read_csv(input_csv)

# Check columns
if 'pdb_id' not in df.columns or 'pdbx_details' not in df.columns:
    raise ValueError("CSV must contain 'pdb_id' and 'pdbx_details' columns.")

# ---------------- Extract reagents for each row ----------------
df['protein'] = ""
df['salt'] = ""
df['buffer'] = ""
df['peg'] = ""

for idx, row in df.iterrows():
    details = row['pdbx_details']
    reagents = extract_reagents(details)
    df.at[idx, 'protein'] = "; ".join(reagents['protein'])
    df.at[idx, 'salt'] = "; ".join(reagents['salt'])
    df.at[idx, 'buffer'] = "; ".join(reagents['buffer'])
    df.at[idx, 'peg'] = "; ".join(reagents['peg'])

# ---------------- Save output ----------------
output_csv = "pdb_results_with_reagents.csv"
df.to_csv(output_csv, index=False)
print(f"Extracted reagents saved to '{output_csv}'")

