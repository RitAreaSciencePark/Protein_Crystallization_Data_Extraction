import os
import json
import pandas as pd
import re

import json

def load_reagents_json(json_path):
    """
    Load the new JSON format with categories and 'name' + 'synonyms',
    and flatten it into a lookup dict:
        reagent_name -> (normalized_name, [classes])
    """
    with open(json_path, "r") as f:
        data = json.load(f)

    lookup = {}

    for category, entries in data.items():
        for entry in entries:
            name = entry["name"].strip()
            # Convert category names to class names (capitalize first letter)
            cls = category.capitalize()

            # Add the main name
            if name not in lookup:
                lookup[name] = (name, [cls])
            else:
                lookup[name][1].append(cls)

            # Add synonyms
            for syn in entry.get("synonyms", []):
                syn = syn.strip()
                if syn not in lookup:
                    lookup[syn] = (name, [cls])
                else:
                    lookup[syn][1].append(cls)

    return lookup


def extract_reagents_with_concentration(text, lookup):
    """
    Detect reagents and extract concentrations (including ranges) from pdbx_details text.
    Handles PEGs and other reagents using a robust lookup.
    
    Returns a dict: {class_name: [list of 'Normalized Reagent Name + conc/range']}
    """
    # 1️⃣ Initialize result dict dynamically from lookup classes
    classes = set(c for c, _ in lookup.values())
    result = {cls: [] for cls in classes}

    if not isinstance(text, str) or not text.strip():
        return result

    # 2️⃣ Normalize text
    text_norm = text.lower()
    text_norm = re.sub(r'[-]', ' ', text_norm)  # replace hyphens with spaces
    text_norm = re.sub(r'\s+', ' ', text_norm)  # collapse multiple spaces

    # 3️⃣ Regex for concentrations
    conc_pattern = re.compile(
        r'(\d+(?:[\.,]\d+)?(?:\s*(?:–|-)\s*\d+(?:[\.,]\d+)?)?\s*(?:%|w/v|m|mm|mol|mmol|mg/ml|g/l|M|mM|uM|µM)?)'
    )

    # 4️⃣ PEG regex (all common formats)
    peg_pattern = re.compile(
        r'(\d+(?:[\.,]\d+)?\s*(?:%|w/v)?\s*peg\s*(?:mme\s*)?\d+[kK]?|\bpeg\d+[kK]?)',
        re.IGNORECASE
    )

    # ---- Process PEGs first ----
    for peg_match in peg_pattern.finditer(text_norm):
        peg_text = peg_match.group(0).strip()
        # Extract concentration if present
        conc_match = conc_pattern.search(peg_text)
        if conc_match:
            conc = conc_match.group(1)
            peg_name_text = peg_text.replace(conc, '').strip()
        else:
            conc = ""
            peg_name_text = peg_text

        # Find normalized name via lookup
        norm_name = None
        cls_name = None
        for syn in lookup.keys():
            if 'peg' in syn and syn in peg_name_text:
                cls_name, norm_name_candidate = lookup[syn]
                norm_name = norm_name_candidate
                break

        # If not found, fallback
        if not norm_name:
            norm_name = peg_name_text.upper().replace(' ', '')
            cls_name = 'precipitants'

        # Ensure class exists in result
        if cls_name not in result:
            result[cls_name] = []

        reagent_str = f"{norm_name} {conc}".strip()
        if reagent_str not in result[cls_name]:
            result[cls_name].append(reagent_str)

    # ---- Process other reagents ----
    # Sort by length to match longest synonym first
    sorted_synonyms = sorted(lookup.keys(), key=lambda x: -len(x))
    for syn in sorted_synonyms:
        if 'peg' in syn:  # skip PEG, already processed
            continue
        cls_name, norm_name = lookup[syn]
        # Ensure class exists
        if cls_name not in result:
            result[cls_name] = []

        for match in re.finditer(r'\b' + re.escape(syn) + r'\b', text_norm):
            start_idx = match.end()
            search_text = text_norm[start_idx:start_idx+30]  # lookahead for concentration
            conc_match = conc_pattern.search(search_text)
            if conc_match:
                reagent_str = f"{norm_name} {conc_match.group(1)}"
            else:
                reagent_str = norm_name
            if reagent_str not in result[cls_name]:
                result[cls_name].append(reagent_str)

    return result

def process_csv(input_csv_path, output_csv_path, json_path):

    lookup = load_reagents_json(json_path)

   # Get all unique classes
    all_classes = set()
    for _, cls_list in lookup.values():
        for cls in cls_list:
            all_classes.add(cls)
    classes = sorted(all_classes)

    df = pd.read_csv(input_csv_path)

    # Create empty columns for each class
    for cls in classes:
        df[cls] = ""

    # Function to extract reagents from pdbx_text
    def extract_reagents_with_concentration(pdbx_text, lookup):
        """
        Dummy parser: extract reagents and assign to classes based on lookup.
        This should be replaced with your actual parsing logic.
        """
        parsed = {cls: [] for cls in classes}
        words = pdbx_text.split(";")  # Example separator
        for word in words:
            word = word.strip()
            if word in lookup:
                norm_name, cls_list = lookup[word]
                for cls in cls_list:
                    parsed[cls].append(norm_name)
        return parsed
  # Process each row
    for idx, row in df.iterrows():
        pdbx_text = row.get("pdbx_details", "")
        parsed = extract_reagents_with_concentration(pdbx_text, lookup)
        for cls in classes:
            df.at[idx, cls] = "; ".join(parsed[cls])

    # Save CSV
    os.makedirs(os.path.dirname(output_csv_path), exist_ok=True)
    df.to_csv(output_csv_path, index=False)
    print(f"✔ Parsed reagents with concentrations saved to {output_csv_path}")