import requests
import gemmi


def fetch_ligands(pdb_id):
    """
    Fetch ligand IDs for a PDB entry.

    1. Try RCSB GraphQL API
    2. If none found, fallback to mmCIF parsing

    Returns:
        list of ligand IDs
    """

    pdb_id = pdb_id.upper()

    # ---------- GraphQL ----------
    url = "https://data.rcsb.org/graphql"

    query = """
    query getLigands($id: String!) {
      entry(entry_id: $id) {
        nonpolymer_entities {
          rcsb_nonpolymer_entity_container_identifiers {
            nonpolymer_comp_id
          }
        }
      }
    }
    """

    try:

        response = requests.post(
            url,
            json={"query": query, "variables": {"id": pdb_id}},
            timeout=15
        )

        response.raise_for_status()

        data = response.json()

        entities = data.get("data", {}).get("entry", {}).get("nonpolymer_entities", [])

        ligands = [
            e["rcsb_nonpolymer_entity_container_identifiers"]["nonpolymer_comp_id"]
            for e in entities
            if e.get("rcsb_nonpolymer_entity_container_identifiers", {}).get("nonpolymer_comp_id")
        ]

        if ligands:
            return list(dict.fromkeys(ligands))

    except Exception as e:
        print(f"GraphQL error: {e}")

    # ---------- mmCIF fallback ----------
    try:

        cif_url = f"https://files.rcsb.org/download/{pdb_id}.cif"

        cif_text = requests.get(cif_url, timeout=15).text

        doc = gemmi.cif.read_string(cif_text)

        block = doc.sole_block()

        comp_ids = block.find_values("_pdbx_entity_instance_feature.comp_id")

        comp_ids = [c for c in comp_ids if c not in ("?", ".", "", None)]

        return list(dict.fromkeys(comp_ids))

    except Exception as e:
        print(f"mmCIF parsing error: {e}")

    return []


# -------- run the function --------
if __name__ == "__main__":

    pdb_id = input("Enter PDB ID: ")

    ligands = fetch_ligands(pdb_id)

    print("Ligands:", ligands)