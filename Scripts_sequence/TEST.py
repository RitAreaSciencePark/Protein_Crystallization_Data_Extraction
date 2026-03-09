import os
import json
import pandas as pd
import re

# -------------------------------
# Classification mapping
# -------------------------------
def map_classification(raw_class):
    raw_class_lower = raw_class.lower()
    if "precipitant" in raw_class_lower or "peg" in raw_class_lower or "mpd" in raw_class_lower or "alcohol" in raw_class_lower:
        return "precipitant"
    elif "buffer" in raw_class_lower or "tris" in raw_class_lower or "hepes" in raw_class_lower:
        return "buffer"
    elif "salt" in raw_class_lower or "chloride" in raw_class_lower or "sulfate" in raw_class_lower or "phosphate" in raw_class_lower or "acetate" in raw_class_lower:
        return "salts"
    else:
        return "additive"

# -------------------------------
# Load JSON dictionary
# -------------------------------
def load_json_classification(json_path):
    with open(json_path, "r") as f:
        data = json.load(f)

    lookup = {}
    for entry in data.get("crystallization_reagents", []):
        canonical_name = entry["canonical"]
        broad_class = map_classification(entry["classification"])
        lookup[canonical_name.lower()] = (canonical_name, broad_class)
        for syn in entry.get("synonyms", []):
            lookup[syn.lower()] = (canonical_name, broad_class)
    return lookup

# -------------------------------
# Extract reagents with concentrations
# -------------------------------
def extract_reagents_classified(text, lookup):
    """
    Extract reagents and concentrations from pdbx_details.
    - PEG: robust extraction (10% (w/v) PEG3350, 35% w/v PEG 8K, etc.)
    - Other reagents: numeric concentration immediately before the compound with units preserved exactly
    - Supports: M, mM, µM, uM, nM, %, mg/ml, mg/l, g/l, mol/l
    - Preserves exact compound text
    """
    result = {}
    used_spans = []

    if not isinstance(text, str):
        return result

    text_norm = text.lower()
    peg_detector = re.compile(r'peg|polyethylene glycol', re.IGNORECASE)

    # PEG patterns
    conc_pattern_peg = re.compile(r'(\d+(?:\.\d+)?\s*%\s*(?:\(\s*(?:w/v|v/v|w/w)\s*\)|w/v|v/v|w/w)?)', re.IGNORECASE)
    conc_pattern_peg_fallback = re.compile(r'(\d+(?:\.\d+)?)', re.IGNORECASE)

    # General concentration pattern for other reagents
    conc_pattern_general = re.compile(r'(\d+(?:\.\d+)?\s*(?:%|M|mM|µM|uM|nM|mg/ml|mg/l|g/l|mol/l))', re.IGNORECASE)

    sorted_synonyms = sorted(lookup.keys(), key=len, reverse=True)

    for syn in sorted_synonyms:
        norm_name, cls_name = lookup[syn]
        pattern = re.compile(r'\b{}\b'.format(re.escape(syn)), re.IGNORECASE)

        for match in pattern.finditer(text_norm):
            span = match.span()
            if any(not (span[1] <= s[0] or span[0] >= s[1]) for s in used_spans):
                continue
            used_spans.append(span)

            original_reagent = text[span[0]:span[1]]
            before = text[max(0, span[0]-50):span[0]]

            if peg_detector.search(original_reagent):
                conc_matches = list(conc_pattern_peg.finditer(before))
                if conc_matches:
                    conc = conc_matches[-1].group(1)
                else:
                    conc_matches = list(conc_pattern_peg_fallback.finditer(before))
                    conc = conc_matches[-1].group(1) if conc_matches else ""
                if conc and "%" not in conc and conc_matches and conc_matches[-1].re == conc_pattern_peg:
                    conc += "%"
            else:
                conc_matches = list(conc_pattern_general.finditer(before))
                conc = conc_matches[-1].group(1) if conc_matches else ""

            reagent = f"{conc + ' ' if conc else ''}{original_reagent}".strip()
            result.setdefault(cls_name, [])
            if reagent.lower() not in [r.lower() for r in result[cls_name]]:
                result[cls_name].append(reagent)

    return result

# -------------------------------
# Process CSV
# -------------------------------
def process_csv_classified(input_csv_path, output_csv_path, json_path):
    lookup = load_json_classification(json_path)
    df = pd.read_csv(input_csv_path)

    classes = ["precipitant", "buffer", "salts", "additive"]
    reagent_dicts = []

    for _, row in df.iterrows():
        text = row.get("pdbx_details", "")
        reagents = extract_reagents_classified(text, lookup)
        reagent_dicts.append(reagents)

    out_df = pd.DataFrame()
    out_df["pdb_id"] = df["PDB_ID"]

    for cls in classes:
        out_df[cls] = ["; ".join(d.get(cls, [])) for d in reagent_dicts]

    os.makedirs(os.path.dirname(output_csv_path), exist_ok=True)
    out_df.to_csv(output_csv_path, index=False)
    print(f"✔ Parsing completed. Output saved to: {output_csv_path}")

# -------------------------------
# Run
# -------------------------------
if __name__ == "__main__":
    json_path = "/home/ruth/Protein_Crystallization_Data_Extraction/Data/reagents.json"
    input_csv_path = "/home/ruth/Protein_Crystallization_Data_Extraction/output/kaiB/kaiB_pdb_mmcif_filtered.csv"
    output_csv_path = "/home/ruth/Protein_Crystallization_Data_Extraction/output/kaiB/kaiB_reagents_parsed.csv"

    process_csv_classified(input_csv_path, output_csv_path, json_path)