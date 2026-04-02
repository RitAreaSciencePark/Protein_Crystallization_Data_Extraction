import requests
from Bio.Align import PairwiseAligner

# --------------------------------------------------
# 1️⃣ Search PDB (get entity_id + RCSB score)
# --------------------------------------------------
def fetch_polymer_entity_ids(sequence):
    url = "https://search.rcsb.org/rcsbsearch/v2/query"

    payload = {
        "query": {
            "type": "terminal",
            "service": "sequence",
            "parameters": {
                "value": sequence,
                "sequence_type": "protein",
                "identity_cutoff": 0.5,
                "evalue_cutoff": 0.1
            }
        },
        "return_type": "polymer_entity",
        "request_options": {
            "paginate": {"start": 0, "rows": 100}
        }
    }

    r = requests.post(url, json=payload, timeout=20)
    r.raise_for_status()

    result_set = r.json().get("result_set", [])

    # return (entity_id, score)
    return [(res["identifier"], res.get("score", 0)) for res in result_set]


# --------------------------------------------------
# 2️⃣ Fetch sequence
# --------------------------------------------------
def fetch_polymer_sequence(entity_id):
    pdb_id, entity_num = entity_id.split("_")
    url = f"https://data.rcsb.org/rest/v1/core/polymer_entity/{pdb_id}/{entity_num}"

    try:
        r = requests.get(url, timeout=10)
        r.raise_for_status()
        data = r.json()

        seq = data.get("entity_poly", {}).get("pdbx_seq_one_letter_code", "")
        return seq.replace("\n", "")

    except Exception:
        return ""


# --------------------------------------------------
# 3️⃣ Compute identity
# --------------------------------------------------
def compute_identity(seq1, seq2):
    aligner = PairwiseAligner()
    score = aligner.score(seq1, seq2)
    return (score / max(len(seq1), len(seq2))) * 100


# --------------------------------------------------
# 4️⃣ Main
# --------------------------------------------------
def search_pdb(sequence):
    print("▶ Searching PDB...")

    hits = fetch_polymer_entity_ids(sequence)
    print(f"✔ Found {len(hits)} hits")

    results = []

    for entity_id, rcsb_score in hits:
        pdb_id = entity_id.split("_")[0]

        seq = fetch_polymer_sequence(entity_id)
        if not seq:
            continue

        identity = compute_identity(sequence, seq)

        results.append((pdb_id, rcsb_score, identity))

    if not results:
        print("❌ No valid results")
        return []

    # Keep best identity per PDB_ID
    best = {}
    for pdb_id, rcsb_score, identity in results:
        if pdb_id not in best or identity > best[pdb_id][1]:
            best[pdb_id] = (rcsb_score, identity)

    # Sort by identity
    sorted_results = sorted(
        [(pdb_id, sc, idt) for pdb_id, (sc, idt) in best.items()],
        key=lambda x: x[2],
        reverse=True
    )

    return sorted_results


# --------------------------------------------------
# ▶ RUN
# --------------------------------------------------
if __name__ == "__main__":
    query_seq = "MSPRKTYILKLYVAGNTPNSVRALKTLKNILEVEFQGVYALKVIDVLKNPQLAEEDKILATPTLAKVLPLPVRRIIGDLSDREKVLIGLDLLYGELQDSDDF"

    results = search_pdb(query_seq)

    print("\nPDB_ID\tRCSB_Score\tIdentity(%)")
    for pdb_id, sc, idt in results[:20]:
        print(f"{pdb_id}\t{sc:.3f}\t\t{idt:.2f}")