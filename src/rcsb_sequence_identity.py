import requests
import pandas as pd

RCSB_SEARCH_URL = "https://search.rcsb.org/rcsbsearch/v2/query"
RCSB_ENTRY_URL = "https://data.rcsb.org/rest/v1/core/entry"


# ─────────────────────────────────────────────
# Detect sequence type
# ─────────────────────────────────────────────
def detect_seq_type(sequence: str):
    seq = sequence.upper()

    if set(seq) <= set("ACDEFGHIKLMNPQRSTVWY"):  # protein letters
        return "protein"
    elif set(seq) <= set("ATCG"):
        return "dna"
    elif set(seq) <= set("AUCG"):
        return "rna"
    else:
        return "protein"


# ─────────────────────────────────────────────
# STEP 1: flexible sequence search
# ─────────────────────────────────────────────
def search_pdb_by_sequence(sequence, max_hits=1000):

    seq_type = detect_seq_type(sequence)

    target_map = {
        "protein": "pdb_protein_sequence",
        "dna": "pdb_dna_sequence",
        "rna": "pdb_rna_sequence"
    }

    target = target_map[seq_type]

    query = {
        "query": {
            "type": "terminal",
            "service": "sequence",
            "parameters": {
                "target": target,
                "value": sequence
            }
        },
        "request_options": {
            "scoring_strategy": "sequence",
            "results_verbosity": "verbose",
            "paginate": {
                "start": 0,
                "rows": max_hits
            }
        },
        "return_type": "polymer_entity"
    }

    r = requests.post(RCSB_SEARCH_URL, json=query)
    r.raise_for_status()
    return r.json().get("result_set", [])


# ─────────────────────────────────────────────
# STEP 2: X-ray filter
# ─────────────────────────────────────────────
def is_xray(pdb_id: str) -> bool:
    try:
        r = requests.get(f"{RCSB_ENTRY_URL}/{pdb_id}", timeout=10)
        if r.status_code != 200:
            return False

        data = r.json()
        method = data.get("exptl", [{}])[0].get("method", "")

        return "X-RAY" in method.upper()

    except Exception:
        return False


# ─────────────────────────────────────────────
# STEP 3: parse hits
# ─────────────────────────────────────────────
def search_rcsb_xray_only(sequence, max_hits=1000):

    raw_hits = search_pdb_by_sequence(sequence, max_hits)

    results = []

    for hit in raw_hits:

        identifier = hit.get("identifier", "")
        score = hit.get("score", None)

        pdb_id = identifier.split("_")[0]
        entity = identifier.split("_")[1] if "_" in identifier else "1"

        seq_identity = None
        evalue = None

        for svc in hit.get("services", []):
            if svc.get("service_type") == "sequence":
                for node in svc.get("nodes", []):
                    for ctx in node.get("match_context", []):
                        seq_identity = ctx.get("sequence_identity")
                        evalue = ctx.get("evalue")

        # ── X-ray filter ──
        if not is_xray(pdb_id):
            continue

        results.append({
            "PDB_ID": pdb_id,
            "Entity": entity,
            "Score": score,
            "Seq_id": seq_identity,
            "E-value": evalue
        })

    return pd.DataFrame(results)


# ─────────────────────────────────────────────
# STEP 4: run + save
# ─────────────────────────────────────────────
def run_and_save (sequence, output_csv_1):

    df = search_rcsb_xray_only(sequence)

    if df.empty:
        print("⚠ No X-ray hits found")
        return df

    df["Score"] = pd.to_numeric(df["Score"], errors="coerce")
    df["Seq_id"] = pd.to_numeric(df["Seq_id"], errors="coerce")

    df = df.sort_values(["Seq_id", "Score"], ascending=False)

    df.to_csv(output_csv_1, index=False)

    print(f"✔ Saved {len(df)} X-ray hits → {output_csv_1}")

    return df