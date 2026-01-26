import requests
import json
import gemmi
import csv
import os
import re

def get_sequence_from_user():
    seq = input(
        "Enter protein sequence (single-letter amino acid code, no FASTA header):\n"
    ).strip().upper()

    # Remove whitespace and line breaks
    seq = re.sub(r"\s+", "", seq)

    if not seq:
        raise ValueError("Sequence cannot be empty")

    # Allowed amino acids for RCSB sequence search
    if not re.fullmatch(r"[ACDEFGHIKLMNPQRSTVWY]+", seq):
        raise ValueError(
            "Sequence contains invalid characters. "
            "Use only standard amino acids (ACDEFGHIKLMNPQRSTVWY)."
        )

    return seq


def run_pdb_search(identity_cutoff=0.9, row_count=20):
    """
    Run a PDB search for a given sequence and extract mmCIF experimental data.
    """
    # Ask for sequence HERE
    seq = get_sequence_from_user()
    # ---------- OUTPUT FILE ----------
    output_csv = "pdb_mmcif_extracted_data.csv"
    ROW_COUNT = 20


#this query is to navigate the PDB database and find the proteins base on the sequences determined by X-ray diffraction ,methods
    query = {
    "query": {
        "type": "group",
        "logical_operator": "and",
        "nodes": [
             {
            "type": "terminal",
            "service": "sequence",
            "parameters": {
                "evalue_cutoff": 1,
                "identity_cutoff": 0.5,
                "sequence_type": "protein",
                "value": seq
            }
            },
            
            
            {
            "type": "terminal",
            "service": "text",
            "parameters": {
                "attribute": "exptl.method",
                "operator": "exact_match",
                "value": "X-RAY DIFFRACTION"
            }
            }
        ]
        },
        "request_options": {
        "paginate": {
            "start": 0,
            "rows": ROW_COUNT
        }
        },
        "return_type": "entry"
    }

    # Load the query string/JSON into a PYTHON dictionary
    #query = json.loads(query)
    query["request_options"]["paginate"]["rows"] = ROW_COUNT

    # Post our query to the system
    result = requests.post("https://search.rcsb.org/rcsbsearch/v2/query",json=query)

    # If you want to know the total number of pages (e.g., to do an automated extraction)
    result.raise_for_status()  # stops if request failed

    # Let's select the first result:

    if result.json()["total_count"] > 0:
        print("First result:", result.json()["result_set"][0])
    else:
        print("No matching structures found.")

    #Get a PDB ID from the search results
    identifier = result.json()["result_set"][0]["identifier"]

    #Fetch detailed metadata for that PDB entry
    entry_data = requests.get(f'https://data.rcsb.org/rest/v1/core/entry/{identifier}').json()

    # Loop over all results and print some information about the experimental conditions
    search_data = result.json()

    if search_data["total_count"] == 0:
        raise ValueError("No matching PDB entries found")

    print("\nExtracting experimental conditions")

    pdb_hits = {hit["identifier"]: hit.get("score") for hit in search_data["result_set"]}

    def extract_pubmed_id(pdb_id):
        """
        Fetch PubMed ID for a given PDB entry.
        """
        entry_data = requests.get(
            f"https://data.rcsb.org/rest/v1/core/entry/{pdb_id}").json()

        pubmed = entry_data.get("rcsb_primary_citation", {}).get("pdbx_database_id_pub_med", "NA")

        return pubmed

    import re

    def extract_protein_concentration(details_text):
        """
        Extract only the protein concentration from pdbx_details.
        Returns the concentration (e.g., '5 mg/mL') or None if not found.
        """
        if not details_text or details_text in {".", "?"}:
            return None

        text = details_text.replace("\n", " ").strip()

        # Match concentrations like 5 mg/mL, 10 mg/ml, 0.5 mg/mL
        pattern = re.compile(r"(\d+(?:\.\d+)?\s*mg/mL)", re.IGNORECASE)

        matches = pattern.findall(text)

        # Return the first match only, which is usually the main protein concentration
        return matches[0] if matches else None

    
    def extract_reagents(details_text):
        """
        Extract reagents (salts, buffers, precipitants, additives) from pdbx_details.
        Excludes protein and protein concentration.
        """
        if not details_text or details_text in {".", "?"}:
            return None

        text = details_text.replace("\n", " ").strip()

        # Remove protein concentration mentions to avoid contamination
        text = re.sub(
            r"\d+(?:\.\d+)?\s*mg/mL\s+[A-Za-z0-9][A-Za-z0-9\-\s]*",
            "",
            text,
            flags=re.IGNORECASE)

        reagent_pattern = re.compile(
            r"""
            (
            \d+(?:\.\d+)?\s*(?:mM|M|%)     # concentration
            \s*
            [A-Za-z][A-Za-z0-9\-\s]*?      # reagent name
            )
            (?=
            , | ; | \bph\b | $             # stop conditions
            )
            """,
            re.IGNORECASE | re.VERBOSE
        )

        reagents = []

        for match in reagent_pattern.findall(text):
            reagent = match.strip()

            # Remove trailing pH info
            reagent = re.sub(r"\s*pH\s*[0-9\.]+", "", reagent, flags=re.IGNORECASE)

            # Reject method words
            if re.search(r"\b(drop|vapor|diffusion|hanging|sitting)\b", reagent, re.I):
                continue

            reagents.append(reagent)

        # Deduplicate while preserving order
        reagents = list(dict.fromkeys(reagents))

        return "; ".join(reagents) if reagents else None


    def parse_from_pdbx_details(details_text):
        """
        Extract crystallization conditions from pdbx_details free text.
        Returns a dictionary with values suitable for CSV output.
        """
        if not details_text or details_text in {".", "?"}:
            return {
                "pH": None,
                "pH_range": None,
                "time": None,
                "seeding": None,
                "buffer": None,
                "peg": None,
                "salt": None,
                "ligand": None,
                "solvent": None
            }

        result = {}
        text = details_text.lower()

        result["protein_concentration"] = extract_protein_concentration(details_text)

        # --- pH range ---
        ph_range_matches = re.findall(r"pH\s*[=:]?\s*([0-9]+(?:\.[0-9]+)?)\s*(?:-|–|to)\s*([0-9]+(?:\.[0-9]+)?)", text)
        pH_range = ", ".join([f"{a}-{b}" for a, b in ph_range_matches]) if ph_range_matches else None

        # --- time ---
        time_matches = re.findall(r"(\d+(?:\.\d+)?)\s*(hours?|hrs?|h|days?|weeks?)", text)
        time = ", ".join(["".join(m) for m in time_matches]) if time_matches else None

        # --- seeding ---
        seeding_matches = re.findall(r"(microseeding|macroseeding|seeding|seeded|unseeded|no seeding)", text)
        seeding = ", ".join(seeding_matches) if seeding_matches else None

        protein_matches = re.findall(r"(\d+(?:\.\d+)?\s*(?:mg/ml|mg\/ml|g/l|µm|um|mm))\s+([a-z0-9\-]+(?:\s+[a-z0-9\-]+)*)", text)
        protein_filtered = []
        for conc, prot in protein_matches:
            if not re.search(r"\b(tris|hepes|mes|mops|cacodylate|acetate|phosphate|citrate|nacl|kcl|mgcl2|cacl2|zncl2|ammonium sulfate|lithium sulfate|sodium sulfate|potassium chloride)\b", prot):
                protein_filtered.append(f"{conc} {prot}")
        protein_concentration = ", ".join(protein_filtered) if protein_filtered else None

        # --- buffer (unique, no repetition) ---
        buffer_matches = re.findall(r"\b(tris|hepes|mes|mops|bis[-\s]?tris|cacodylate|acetate|phosphate|citrate)\b\s*(hcl)?", text)
        buffer = ", ".join(sorted(set(" ".join(filter(None, m)).strip() for m in buffer_matches))) if buffer_matches else None

        # --- PEG (with % and unique) ---
        peg_matches = re.findall(r"(\d{1,3}\s*%\s*.*?peg\s*\d{3,4}|\d{1,3}\s*%\s*.*?polyethylene glycol\s*\d{3,4})", text)
        peg = ", ".join(sorted(set([m.replace("(w/v)", "").strip() for m in peg_matches]))) if peg_matches else None

        # --- salt (unique, no repetition) ---
        salt_matches = re.findall(r"\b(licl|nacl|kcl|mgcl2|cacl2|zncl2|ammonium sulfate|lithium sulfate|sodium sulfate|potassium chloride)\b", text)
        salt = ", ".join(sorted(set(salt_matches))) if salt_matches else None

        # --- ligand ---
        ligand_matches = re.findall(r"([a-z0-9\-]+)\s*(?:ligand|inhibitor|substrate|cofactor)", text)
        ligand = ", ".join(sorted(set(ligand_matches))) if ligand_matches else None

        # --- solvent ---
        solvent_matches = re.findall(r"\b(water|h2o|ethanol|isopropanol|mpd|glycerol|dmso)\b", text)
        solvent = ", ".join(sorted(set(solvent_matches))) if solvent_matches else None


        return {
            "pH_range": pH_range,
            "time": time,
            "seeding": seeding,
            "protein_concentration": protein_concentration,
            "buffer": buffer,
            "peg": peg,
            "salt": salt,
            "ligand": ligand,
            "solvent": solvent
        }
    def extract_ph_from_details(details_text):
        """
        Extract a single pH value from pdbx_details text.
        """
        if not details_text or details_text in {".", "?"}:
            return None

        text = details_text.lower()

        match = re.search(
            r"\bph\s*[=:]?\s*([0-9]+(?:\.[0-9]+)?)",
            text
        )

        return match.group(1) if match else None



    def extract_mmcif_info(pdb_id, score):
        pdb_id = pdb_id.upper()
        pubmed = extract_pubmed_id(pdb_id)

        # ---- Fetch mmCIF file directly from PDB ----
        url = f"https://files.rcsb.org/view/{pdb_id}.cif"
        response = requests.get(url)
        response.raise_for_status()  # Raise error if download fails

        # ---- Parse mmCIF directly from text ----
        doc = gemmi.cif.read_string(response.text)
        block = doc.sole_block()  # Get the main data block

        pdbx_details = block.find_value("_exptl_crystal_grow.pdbx_details")
        pdbx_ph_range = block.find_value("_exptl_crystal_grow.pdbx_pH_range")
        ph_structured = block.find_value("_exptl_crystal_grow.pH")

        extracted = parse_from_pdbx_details(pdbx_details)

        # Prefer structured mmCIF, fallback to free text
        if pdbx_ph_range in {None, ".", "?"}:
            pdbx_ph_range = extracted["pH_range"]

        # --- Fallback extraction ---
        if ph_structured in {None, ".", "?"}:
            ph_final = extract_ph_from_details(pdbx_details)
        else:
            ph_final = ph_structured
        
        protein_concentration = extract_protein_concentration(pdbx_details)
        reagents = extract_reagents(pdbx_details)

        info = {
            "PDB_ID": block.find_value("_entry.id"),
            "score": score,
            "Resolution": block.find_value("_refine.ls_d_res_high"),
            "pubmed_id": pubmed,
            "method": block.find_value("_exptl_crystal_grow.method"),
            "temp": block.find_value("_exptl_crystal_grow.temp"),
            "pH": ph_final,
            "seeding": extracted["seeding"],
            "time": extracted["time"],
            "pdbx_pH_range": pdbx_ph_range,
            "pH_range_text": extracted["pH_range"],
            "protein_concentration": protein_concentration,
            "buffer": extracted["buffer"],
            "peg": extracted["peg"],
            "salt": extracted["salt"],
            "ligand": extracted["ligand"],
            "solvent": extracted["solvent"],
            "reagents": reagents
            }
        
        return info


    output_csv = "pdbx_details_extracted.csv"

    with open(output_csv, "w", newline="") as f:
        writer = csv.DictWriter(
            f,
            fieldnames=["PDB_ID", "score", "Resolution", "pubmed_id", 
            "method", "pH", "seeding", "temp", "time", "pdbx_pH_range", "pH_range_text", "protein_concentration", "reagents", "buffer", "peg", "salt", "ligand", "solvent"])
        writer.writeheader()

        for pdb_id, score in pdb_hits.items():
            try:
                row= extract_mmcif_info(pdb_id, score)
                writer.writerow(row)
            except Exception as e:
                print(f"Failed for {pdb_id}: {e}")

    print(f"\nCSV file saved as: {os.path.abspath(output_csv)}")

if __name__ == "__main__":
    run_pdb_search()
