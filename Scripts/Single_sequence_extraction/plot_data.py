import pandas as pd
import matplotlib.pyplot as plt
import matplotlib.lines as mlines
import numpy as np
import re

# ---------- Load data ----------
df = pd.read_csv("pdbx_details_extracted.csv")

# ---------- Clean numeric ----------
df["pH"] = pd.to_numeric(df["pH"], errors="coerce")
df["temp"] = pd.to_numeric(df["temp"], errors="coerce")

df = df.dropna(subset=["pH", "temp", "method", "PDB_ID"])

# ---------- Keep hanging / sitting ----------
df = df[df["method"].str.contains("hanging|sitting", case=False, na=False)]

# ---------- Parse pH range ----------
def parse_ph_range(ph_range):
    if pd.isna(ph_range):
        return None, None
    m = re.match(r"([0-9.]+)\s*-\s*([0-9.]+)", str(ph_range))
    if m:
        return float(m.group(1)), float(m.group(2))
    return None, None

df[["pH_low", "pH_high"]] = df["pdbx_pH_range"].apply(
    lambda x: pd.Series(parse_ph_range(x))
)

df["pH_err_low"] = df["pH"] - df["pH_low"]
df["pH_err_high"] = df["pH_high"] - df["pH"]

# ---------- Marker mapping ----------
def get_marker(method):
    if "hanging" in method.lower():
        return "^"
    elif "sitting" in method.lower():
        return "s"
    return "o"

df["marker"] = df["method"].apply(get_marker)

# ---------- Assign colors per PDB_ID ----------
pdb_ids = sorted(df["PDB_ID"].unique())
colors = plt.cm.tab20(np.linspace(0, 1, len(pdb_ids)))
color_map = dict(zip(pdb_ids, colors))

# ---------- Plot ----------
plt.figure(figsize=(10, 7))

for _, row in df.iterrows():
    plt.errorbar(
        row["temp"],
        row["pH"],
        yerr=[[row["pH_err_low"]], [row["pH_err_high"]]]
        if pd.notna(row["pH_err_low"]) and pd.notna(row["pH_err_high"])
        else None,
        fmt=row["marker"],
        color=color_map[row["PDB_ID"]],
        ecolor=color_map[row["PDB_ID"]],
        markersize=8,
        markeredgecolor="black",
        capsize=3,
        elinewidth=1
    )

# ---------- Axes ----------
plt.xlabel("Temperature (K)")
plt.ylabel("pH")
plt.title("pH vs Temperature (K)\nshape = method | color = PDB ID")
plt.grid(True)

# ---------- X-axis: show temperature values of points ----------
unique_temps = sorted(df["temp"].unique())
plt.xticks(unique_temps, [f"{t:.1f}" for t in unique_temps], rotation=45)

# ---------- Legends ----------
pdb_legend = [
    mlines.Line2D([], [], color=color_map[pdb_id], marker='o',
                  linestyle='None', markersize=8, label=pdb_id)
    for pdb_id in pdb_ids
]

method_legend = [
    mlines.Line2D([], [], color='black', marker='^',
                  linestyle='None', markersize=8, label='Hanging drop'),
    mlines.Line2D([], [], color='black', marker='s',
                  linestyle='None', markersize=8, label='Sitting drop')
]

legend1 = plt.legend(handles=pdb_legend, title="PDB ID (color)",
                     bbox_to_anchor=(1.02, 1), loc="upper left")
plt.gca().add_artist(legend1)

plt.legend(handles=method_legend, title="Method (shape)",
           bbox_to_anchor=(1.02, 0.45), loc="upper left")

# ---------- Save as PDF ----------
output_pdf = "pH_vs_Temperature_K_xticks.pdf"
plt.tight_layout()
plt.savefig(output_pdf)
plt.show()

print(f"Plot saved as: {output_pdf}")
