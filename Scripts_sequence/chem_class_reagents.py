import os
import json
import pandas as pd
import re
import json

def load_json(json_path):
    """Load a simple synonym -> normalized name JSON."""
    with open(json_path, "r") as f:
        lookup = json.load(f)

    # Make all keys lowercase for matching
    lookup = {k.lower(): v for k, v in lookup.items()}
    return lookup

def extract_reagents_with_concentration_simple(text, lookup):
    """
    Extract reagents and concentrations from pdbx_details text.
    Returns a list of strings: 'Normalized Name + concentration'
    """
    reagents = []
    if not text:
        return reagents

    text_norm = text.lower()

    # Regex for concentration: number + optional decimal + optional unit
    conc_pattern = re.compile(r'(\d+(?:[\.,]\d+)?\s*(?:%|M|mM|uM|µM|mg/ml|g/l)?)')

    # Sort lookup keys by length (longest first) to match first
    sorted_keys = sorted(lookup.keys(), key=lambda x: -len(x))

    for syn in sorted_keys:
        # Match whole word
        if re.search(r'\b{}\b'.format(re.escape(syn)), text_norm):
            # Look for concentration within ~20 chars after match
            match_pos = re.search(r'\b{}\b'.format(re.escape(syn)), text_norm).end()
            nearby_text = text_norm[match_pos:match_pos+20]
            conc_match = conc_pattern.search(nearby_text)
            conc_str = conc_match.group(1) if conc_match else ''
            reagent_str = lookup[syn]
            if conc_str:
                reagent_str += f" {conc_str}"
            if reagent_str not in reagents:
                reagents.append(reagent_str)
    return reagents

def process_csv(input_csv_path, output_csv_path, json_path):
    lookup = load_json(json_path)
    df = pd.read_csv(input_csv_path)

    # Add a column for reagents with concentrations
    df['reagents_parsed'] = ""

    for idx, row in df.iterrows():
        pdbx_text = row.get("pdbx_details", "")
        reagents_list = extract_reagents_with_concentration_simple(pdbx_text, lookup)
        df.at[idx, 'reagents_parsed'] = "; ".join(reagents_list)

    os.makedirs(os.path.dirname(output_csv_path), exist_ok=True)
    df.to_csv(output_csv_path, index=False)
    print(f"✔ Parsed reagents saved to {output_csv_path}")