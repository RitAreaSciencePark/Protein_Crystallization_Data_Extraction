import pandas as pd

# ===== 1. File paths =====
input_file = r"/u/mdmc/rnananja/Ruthprojet1/Protein_Crystallization_Data_Extraction/output/kaiB/kaiB_pdb_mmcif_filtered.csv"
output_file = r"/u/mdmc/rnananja/Ruthprojet1/Protein_Crystallization_Data_Extraction/output/kaiB/crystallization_data.xlsx"

# ===== 2. Load CSV =====
df = pd.read_csv(input_file)

# ===== 3. Columns to group =====
cocktail_cols = ["pH", "COMPOUNDS", "temp", "method", "ligands"]

# ===== 4. Create hierarchical columns =====
new_columns = []
for col in df.columns:
    if col in cocktail_cols:
        new_columns.append(("Crystallization_cocktails", col))
    else:
        new_columns.append((col, ""))

df.columns = pd.MultiIndex.from_tuples(new_columns)

# ===== 5. Optional: move grouped columns to the end =====
cocktail = df["Crystallization_cocktails"]
others = df.drop(columns="Crystallization_cocktails", level=0)
df = pd.concat([others, cocktail], axis=1)

# ===== 6. Save Excel file =====
df.to_excel(output_file, index=False)

print("File successfully saved to:")
print(output_file)