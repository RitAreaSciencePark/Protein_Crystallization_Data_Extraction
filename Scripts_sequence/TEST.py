import requests

def check_ligands(pdb_id):
    url = f"https://www.rcsb.org/ligand-validation/{pdb_id}"
    try:
        response = requests.get(url, timeout=5)
        # The page always returns 200, so check for the "Coerced Null value" message
        if "Coerced Null value" in response.text:
            ligands = "none"
        else:
            ligands = f"yes ({url})"
    except requests.RequestException:
        # Any connection issue, timeout, etc.
        ligands = "none"
    return ligands

# Example usage
pdb_id = "4KSO" 
ligand_status = check_ligands(pdb_id)
print(f"PDB ID {pdb_id} ligands: {ligand_status}")