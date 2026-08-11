"""
plot.py
=======

Generates the crystallization-condition plots (pH vs Temperature, pH vs
PEG concentration) and the colored PDF summary table for one protein's
search results.

Reads
-----
`Output_compounds.csv` from the protein's output folder (the file
compound_extraction.process_csv() produces), with columns:

    pdb_id, entity_id, sequence_identity, evalue, score, Resolution,
    Pubmed_id, Polymer, Assembly, Method, pH, Temp, pdbx_details,
    pdbx_pH_range, Non_polymers, compound

`compound` holds the aggregated string compound_extraction.py produces,
e.g. "'PEG 8000' (10% w/v), 'Tris-HCl' (0.1 M)" -- there's no separate
PEG_con column anymore, so PEG concentration for the pH-vs-PEG plot is
parsed straight out of that string.

`Non_polymers` holds an RCSB 3D-view URL (or nothing) rather than ligand
names -- it's shown in the table under "3D View".

Writes (into the same protein output folder)
----------------------------------------------
    {protein_name}_TEMP.png
    {protein_name}_PEG.png
    {protein_name}_Cryst_cocktail_Table.pdf
    Grouped_conditions.csv   <- one row per unique condition, with a
                                `row_id` column shared between the plot
                                point labels and the PDF table rows.

On clickability
----------------
Matplotlib PNG/PDF output is static -- there's no click-event runtime in
an image or PDF, so true "click a point, jump to its table row" has to
live in the interactive web app (Plotly/Django), not here. What this
script *can* do, and does, is make every plotted point and every table
row share the same `row_id`, so the web app's data loader (see
crystal_explorer/viewer/data.py) can read Grouped_conditions.csv and
wire up that interactivity on top of the same grouping logic used here.
"""

import os
import re
import textwrap

import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import matplotlib.lines as mlines
import matplotlib.patheffects as pe


# --------------------------------------------------------------------------- #
# Parsing the aggregated `compound` string
# --------------------------------------------------------------------------- #
# Matches: 'compound name' (concentration) -- concentration may be "nan"
_COMPOUND_ENTRY_RE = re.compile(r"'([^']+)'\s*\(([^)]*)\)")


def parse_compound_string(compound_str):
    """Parse compound_extraction.py's aggregated string into a list of
    (name, concentration) tuples. Returns [] for NaN/empty input."""
    if pd.isna(compound_str) or not isinstance(compound_str, str):
        return []
    return _COMPOUND_ENTRY_RE.findall(compound_str)


def extract_peg_percent(compound_str):
    """Pull the first PEG entry's percentage concentration out of the
    aggregated compound string, e.g. "'PEG 8000' (10% w/v)" -> 10.0.
    Returns np.nan if there's no PEG entry or it has no numeric % value."""
    for name, concentration in parse_compound_string(compound_str):
        if name.strip().upper().startswith("PEG"):
            m = re.search(r"(\d+(?:\.\d+)?)\s*%", concentration)
            if m:
                return float(m.group(1))
    return np.nan


# --------------------------------------------------------------------------- #
# Column widths helper (unchanged)
# --------------------------------------------------------------------------- #
def compute_col_widths(df, scale=0.01):
    """Compute column widths for a matplotlib table based on text length."""
    widths = {}
    for col in df.columns:
        max_len = max(df[col].astype(str).apply(len).max(), len(col))
        widths[col] = max_len * scale
    return widths


# --------------------------------------------------------------------------- #
# MAIN FUNCTION
# --------------------------------------------------------------------------- #
def run_plot(protein_dir, protein_name, compounds_csv_name="Output_compounds.csv"):
    """
    Read {protein_dir}/{compounds_csv_name} and generate the TEMP plot,
    PEG plot, and colored PDF table, all saved into protein_dir.
    """
    output_csv_file = os.path.join(protein_dir, compounds_csv_name)
    print(f"Loading CSV: {output_csv_file}")
    df = pd.read_csv(output_csv_file)

    # ---- Normalize column names to what the rest of this function expects ----
    # (pdb_sequence_search.py / compound_extraction.py use lowercase/different
    # names than this script's original column assumptions.)
    rename_map = {
        "pdb_id": "PDB_ID",
        "score": "Score",
        "sequence_identity": "Seq_id",
    }
    df = df.rename(columns={k: v for k, v in rename_map.items() if k in df.columns})

    for required in ["PDB_ID", "Score", "Seq_id", "Pubmed_id", "Polymer", "Assembly",
                      "Method", "pH", "Temp", "pdbx_pH_range", "compound"]:
        if required not in df.columns:
            df[required] = np.nan

    # "Ligands" is now an RCSB 3D-view URL (Non_polymers), not a ligand list.
    df["Ligands"] = df["Non_polymers"] if "Non_polymers" in df.columns else np.nan

    # --------------------------------------------------------
    # Clean method column
    # --------------------------------------------------------
    df["Method"] = (
        df["Method"]
        .fillna("unspecified")
        .astype(str)
        .str.lower()
        .str.replace(r"[^\w\s]", "", regex=True)
        .str.replace(r"\s+", " ", regex=True)
        .str.strip()
    )
    df.loc[df["Method"] == "", "Method"] = "unspecified"

    # --------------------------------------------------------
    # Numeric conversion
    # --------------------------------------------------------
    df["pH"] = pd.to_numeric(df["pH"], errors="coerce")
    df["Temp"] = pd.to_numeric(df["Temp"], errors="coerce")
    df["Score"] = pd.to_numeric(df["Score"], errors="coerce")
    df["Score"] = df["Score"].clip(0.5, 1)

    # --------------------------------------------------------
    # Parse pH range
    # --------------------------------------------------------
    def parse_ph_range(ph_range):
        if pd.isna(ph_range):
            return (np.nan, np.nan)
        m = re.search(r"(\d+(?:\.\d+)?)\s*[-\u2013]\s*(\d+(?:\.\d+)?)", str(ph_range))
        if m:
            return (float(m.group(1)), float(m.group(2)))
        return (np.nan, np.nan)

    ranges = df["pdbx_pH_range"].apply(parse_ph_range)
    df["pH_low"] = ranges.apply(lambda x: x[0])
    df["pH_high"] = ranges.apply(lambda x: x[1])

    # --------------------------------------------------------
    # Compute plotting pH
    # --------------------------------------------------------
    def compute_plot_ph(row):
        pH = row["pH"] if pd.notna(row["pH"]) else None
        ph_low = row["pH_low"] if pd.notna(row["pH_low"]) else None
        ph_high = row["pH_high"] if pd.notna(row["pH_high"]) else None

        if pH is not None:
            return (float(pH), 0.0, 0.0)
        if ph_low is not None and ph_high is not None:
            ph_low, ph_high = float(ph_low), float(ph_high)
            ph_mid = (ph_low + ph_high) / 2
            return (ph_mid, ph_mid - ph_low, ph_high - ph_mid)
        return (np.nan, 0.0, 0.0)

    result = df.apply(compute_plot_ph, axis=1)
    result = pd.DataFrame(result.tolist(), columns=["plot_pH", "err_low", "err_high"])
    df = pd.concat([df, result], axis=1)

    df["has_ph"] = ~df["plot_pH"].isna()
    df["has_temp"] = ~df["Temp"].isna()

    valid_ph = df["plot_pH"].dropna()
    valid_temp = df["Temp"].dropna()

    ph_min = 4.5
    ph_max = valid_ph.max() if not valid_ph.empty else ph_min + 1
    temp_min = 275
    temp_max = valid_temp.max() if not valid_temp.empty else temp_min + 10

    no_ph_y = ph_min - 0.5
    no_temp_x = temp_min - 4

    df["plot_pH_numeric"] = df["plot_pH"]
    df["pH"] = df["plot_pH_numeric"]

    df["pH_plot"] = df["plot_pH_numeric"].copy()
    df.loc[df["plot_pH_numeric"].isna() | (df["plot_pH_numeric"] < ph_min), "pH_plot"] = no_ph_y

    df["temp_plot"] = df["Temp"].copy()
    df.loc[df["Temp"].isna() | (df["Temp"] < temp_min), "temp_plot"] = no_temp_x

    temp_min_tick = int(np.floor(temp_min / 5) * 5)
    temp_max_tick = 0 if pd.isna(temp_max) else int(np.ceil(temp_max / 5) * 5)

    ph_min_tick = np.floor(ph_min * 2) / 2
    ph_max_tick = np.ceil(ph_max * 2) / 2

    # --------------------------------------------------------
    # PEG concentration, parsed from the aggregated `compound` string
    # --------------------------------------------------------
    df["PEG_con_plot"] = df["compound"].apply(extract_peg_percent)

    # --------------------------------------------------------
    # Marker assignment
    # --------------------------------------------------------
    method_marker_map = {
        "vapor diffusion hanging drop": "s",
        "hanging drop": "s",
        "hanging drop vapor diffusion": "s",
        "vapor diffusion sitting drop": "^",
        "sitting drop": "^",
        "sitting drop vapor diffusion": "^",
        "vapor diffusion": "o",
        "unspecified": "X",
        "": "X",
    }
    fallback_markers = ["D", "v", "P", "*", "<", ">"]
    used_fallback = {}

    def assign_marker(method):
        method = str(method).lower().strip()
        if method in method_marker_map:
            return method_marker_map[method]
        if method not in used_fallback:
            used_fallback[method] = len(used_fallback) % len(fallback_markers)
        return fallback_markers[used_fallback[method]]

    df["marker"] = df["Method"].apply(assign_marker)

    cmap = plt.cm.viridis
    norm = plt.Normalize(0.5, 1)

    # --------------------------------------------------------
    # Group into unique conditions, merging PDB IDs, and assign row_id
    # --------------------------------------------------------
    def all_conditions_with_merged_pdb(df):
        """One row per unique (Pubmed_id, Method, pH, compound) condition,
        with PDB IDs sharing that condition merged into one string, and a
        `row_id` assigned -- shared by the plot point label and the PDF
        table row, and intended to also be shared with the web app."""
        cond_cols = ["Pubmed_id", "Method", "plot_pH_numeric", "compound"]

        key_cols = df[cond_cols].fillna("").astype(str)
        condition_key = key_cols.agg("||".join, axis=1)

        df_with_key = df.copy()
        df_with_key["_condition_key"] = condition_key

        merged_rows = []
        for _, group in df_with_key.groupby("_condition_key"):
            merged_pdb = ", ".join(group["PDB_ID"].astype(str).tolist())
            row = group.iloc[0].copy()
            row["PDB_ID"] = merged_pdb
            merged_rows.append(row)

        merged_df = pd.DataFrame(merged_rows)
        if merged_df.empty:
            return merged_df

        for col in ["Score", "Seq_id", "pH", "Temp", "plot_pH_numeric", "PEG_con_plot"]:
            if col in merged_df.columns:
                merged_df[col] = pd.to_numeric(merged_df[col], errors="coerce")

        merged_df = merged_df.sort_values("Score", ascending=False).reset_index(drop=True)
        merged_df.drop(columns=["_condition_key"], inplace=True)
        merged_df["row_id"] = merged_df.index
        return merged_df

    grouped = all_conditions_with_merged_pdb(df)

    # Persist the grouped/row_id table for the web app to consume later.
    grouped_csv_path = os.path.join(protein_dir, "Grouped_conditions.csv")
    grouped.to_csv(grouped_csv_path, index=False)
    print(f"Grouped conditions (with row_id) saved to: {grouped_csv_path}")

    # ========================================================
    # 1. FULL PLOT pH vs Temp (with Score coloring)
    # ========================================================
    fig, ax_temp = plt.subplots(figsize=(10, 6))

    for _, row in grouped.iterrows():
        x = row["temp_plot"] if pd.notna(row.get("temp_plot")) else no_temp_x
        y = row["pH_plot"] if pd.notna(row.get("pH_plot")) else no_ph_y

        ax_temp.errorbar(
            x, y,
            yerr=[[row["err_low"]], [row["err_high"]]] if row["has_ph"] else None,
            fmt=row["marker"],
            color=cmap(norm(row["Score"] if pd.notna(row["Score"]) else 0.5)),
            ecolor=cmap(norm(row["Score"] if pd.notna(row["Score"]) else 0.5)),
            markeredgecolor="black",
            markersize=9,
            capsize=3,
        )
        # row_id is embedded in the point label so a point can be traced
        # back to its Grouped_conditions.csv / PDF table row by eye or by
        # a downstream tool, even from the static PNG.
        ax_temp.text(x, y, f"{row['PDB_ID']} [#{row['row_id']}]", fontsize=7,
                     path_effects=[pe.withStroke(linewidth=2, foreground="white")])

    ax_temp.axhline(no_ph_y + 0.25, linestyle="--", color="black")
    ax_temp.axvline(no_temp_x + 2, linestyle="--", color="black")

    ax_temp.set_xlim(no_temp_x - 2, valid_temp.max() + 5 if not valid_temp.empty else no_temp_x + 20)
    xticks = [no_temp_x] + list(range(temp_min_tick, temp_max_tick + 5, 10))
    ax_temp.set_xticks(xticks)
    ax_temp.set_xticklabels(["No Temp"] + [str(t) for t in xticks[1:]])

    ax_temp.set_ylim(no_ph_y - 0.25, (valid_ph.max() if not valid_ph.empty else ph_min) + 0.25)
    yticks = [no_ph_y] + list(np.arange(ph_min_tick, ph_max_tick + 0.5, 1))
    ax_temp.set_yticks(yticks)
    ax_temp.set_yticklabels(["No pH"] + [f"{t:.1f}" for t in yticks[1:]])

    ax_temp.set_xlabel("Temperature (K)")
    ax_temp.set_ylabel("pH")
    ax_temp.set_title(f"pH vs Temperature (K)\n {protein_name}", fontsize=14, fontweight="bold")

    unique_methods = grouped["Method"].unique()
    legend_items = []
    used_markers = {}
    for method in unique_methods:
        marker = assign_marker(method)
        if marker not in used_markers:
            legend_items.append(
                mlines.Line2D([], [], color="black", marker=marker, linestyle="None",
                               markersize=8, label=str(method).title()))
            used_markers[marker] = True

    ax_temp.legend(handles=legend_items, title="Method", loc="lower center",
                    bbox_to_anchor=(0.5, -0.25), ncol=3, frameon=False)

    sm = plt.cm.ScalarMappable(cmap=cmap, norm=norm)
    sm.set_array([])
    cbar = fig.colorbar(sm, ax=ax_temp, pad=0.02)
    cbar.set_label("Score", rotation=270, labelpad=15)
    fixed_ticks = [0.5, 0.6, 0.7, 0.8, 0.9, 1.0]
    cbar.set_ticks(fixed_ticks)
    cbar.set_ticklabels([f"{t:.1f}" for t in fixed_ticks])

    temp_png = os.path.join(protein_dir, f"{protein_name}_TEMP.png")
    fig.savefig(temp_png, dpi=300, bbox_inches="tight")
    plt.close(fig)

    # ========================================================
    # 2. FULL PLOT pH vs PEG concentration (with Score coloring)
    # ========================================================
    fig, ax_peg = plt.subplots(figsize=(10, 6))

    no_peg_x = -5
    valid_peg = grouped["PEG_con_plot"].dropna()

    if valid_peg.empty:
        peg_min_tick, peg_max_tick, peg_max = 0, 10, 10
    else:
        peg_min_tick = int(np.floor(valid_peg.min() / 5) * 5)
        peg_max_tick = int(np.ceil(valid_peg.max() / 5) * 5)
        peg_max = valid_peg.max()

    ax_peg.set_xlim(no_peg_x - 2, peg_max + 5)
    xticks = [no_peg_x] + list(range(peg_min_tick, peg_max_tick + 5, 10))
    ax_peg.set_xticks(xticks)
    ax_peg.set_xticklabels(["No PEG"] + [str(t) for t in xticks[1:]])

    for _, row in grouped.iterrows():
        x = row["PEG_con_plot"] if pd.notna(row["PEG_con_plot"]) else no_peg_x
        y = row["pH_plot"] if pd.notna(row.get("pH_plot")) else no_ph_y

        ax_peg.errorbar(
            x, y,
            yerr=[[row["err_low"]], [row["err_high"]]] if row["has_ph"] else None,
            fmt=row["marker"],
            color=cmap(norm(row["Score"] if pd.notna(row["Score"]) else 0.5)),
            ecolor=cmap(norm(row["Score"] if pd.notna(row["Score"]) else 0.5)),
            markeredgecolor="black",
            markersize=9,
            capsize=3,
        )
        ax_peg.text(x, y, f"{row['PDB_ID']} [#{row['row_id']}]", fontsize=7,
                    path_effects=[pe.withStroke(linewidth=2, foreground="white")])

    ax_peg.axhline(no_ph_y + 0.25, linestyle="--", color="black")
    ax_peg.axvline(no_peg_x + 2, linestyle="--", color="black")
    ax_peg.set_xlim(no_peg_x - 2, peg_max + 5)
    ax_peg.set_xticks(xticks)
    ax_peg.set_xticklabels(["No PEG"] + [str(t) for t in xticks[1:]])

    ax_peg.set_ylim(no_ph_y - 0.25, (valid_ph.max() if not valid_ph.empty else ph_min) + 0.25)
    yticks = [no_ph_y] + list(np.arange(ph_min_tick, ph_max_tick + 0.5, 1))
    ax_peg.set_yticks(yticks)
    ax_peg.set_yticklabels(["No pH"] + [f"{t:.1f}" for t in yticks[1:]])

    ax_peg.set_xlabel("PEG concentration (%)")
    ax_peg.set_ylabel("pH")
    ax_peg.set_title(f"pH vs PEG concentration (%)\n {protein_name}", fontsize=18, fontweight="bold")

    legend_items = []
    used_markers = {}
    for method in unique_methods:
        marker = assign_marker(method)
        if marker not in used_markers:
            legend_items.append(
                mlines.Line2D([], [], color="black", marker=marker, linestyle="None",
                               markersize=10, label=str(method).title()))
            used_markers[marker] = True

    ax_peg.legend(handles=legend_items, title="Method", loc="lower center",
                   bbox_to_anchor=(0.5, -0.25), ncol=3, frameon=False)

    sm = plt.cm.ScalarMappable(cmap=cmap, norm=norm)
    sm.set_array([])
    cbar = fig.colorbar(sm, ax=ax_peg, pad=0.02)
    cbar.set_label("Score", rotation=270, labelpad=15)
    cbar.set_ticks(fixed_ticks)
    cbar.set_ticklabels([f"{t:.1f}" for t in fixed_ticks])

    peg_png = os.path.join(protein_dir, f"{protein_name}_PEG.png")
    fig.savefig(peg_png, dpi=300, bbox_inches="tight")
    plt.close(fig)

    # ========================================================
    # Colored PDF TABLE (row_id shown as the first column)
    # ========================================================
    fig, ax = plt.subplots(figsize=(18, 10))
    ax.axis("off")

    def _clean_cell(value):
        """Blank out NaN/None rather than rendering the literal string 'nan'."""
        return "" if pd.isna(value) else value

    def _clean_pubmed(value):
        if pd.isna(value) or str(value) in ("NA", ""):
            return ""
        if isinstance(value, float) and value.is_integer():
            return str(int(value))
        return str(value)

    table_rows = []
    for _, group in grouped.iterrows():
        table_rows.append([
            str(group["row_id"]),
            group["PDB_ID"],
            f"{group['Score']:.3f}" if pd.notna(group["Score"]) else "",
            f"{group['Seq_id']:.3f}" if pd.notna(group["Seq_id"]) else "",
            _clean_pubmed(group["Pubmed_id"]),
            _clean_cell(group.get("Polymer")),
            _clean_cell(group.get("Assembly")),
            group["Method"],
            _clean_cell(group.get("Ligands")),
            group["plot_pH_numeric"] if pd.notna(group["plot_pH_numeric"]) else "",
            group["Temp"] if pd.notna(group["Temp"]) else "",
            _clean_cell(group.get("compound")),
        ])

    main_header = ["", "", "", "", "", "", "", "", "", "", "", "CRYSTALLIZATION COCKTAILS"]
    sub_headers = ["#", "PDB_ID", "Score", "Seq_id", "Pubmed_id", "Polymer", "Assembly",
                   "Method", "3D View", "pH", "Temp", "Compounds"]
    col_widths = [0.03, 0.14, 0.05, 0.05, 0.08, 0.07, 0.09, 0.10, 0.09, 0.04, 0.04, 0.32]

    tbl = ax.table(
        cellText=[main_header, sub_headers] + table_rows,
        colWidths=col_widths,
        cellLoc="center",
        loc="center",
    )

    for c in range(len(main_header)):
        txt = tbl[0, c].get_text()
        txt.set_fontsize(18)
        txt.set_weight("bold")

    for c in range(len(sub_headers)):
        txt = tbl[1, c].get_text()
        txt.set_weight("bold")
        txt.set_fontsize(12)

    col_colors = {
        "#": "#DDDDDD",
        "PDB_ID": "#98ACBA",
        "Score": "#B7EFBC",
        "Seq_id": "#F17DFB",
        "Pubmed_id": "#D8BC8E",
        "Polymer": "#4DAEC8",
        "Assembly": "#CA94D2",
        "Method": "#9AD7DF",
        "3D View": "#E9A6BC",
        "pH": "#818FD8",
        "Temp": "#CBC8B1",
        "Compounds": "#A7F054",
    }

    for col, col_name in enumerate(sub_headers):
        color = col_colors.get(col_name, "white")
        for row in range(1, len(table_rows) + 2):
            tbl[row, col].set_facecolor(color)

    exp_start = sub_headers.index("pH")
    exp_end = sub_headers.index("Compounds")
    group_color = "#F2EEED"
    for c in range(exp_start, exp_end + 1):
        cell = tbl[0, c]
        cell.set_facecolor(group_color)
        if c == exp_start:
            cell.visible_edges = "LTB"
        elif c == exp_end:
            cell.visible_edges = "RTB"
        else:
            cell.visible_edges = "TB"
        cell.set_linewidth(1.2)

    pdb_start = sub_headers.index("PDB_ID")
    pdb_end = sub_headers.index("3D View")
    for c in range(pdb_start, pdb_end + 1):
        cell = tbl[0, c]
        cell.visible_edges = "RB" if c == pdb_end else "B"
        cell.set_linewidth(1.2)

    wrap_widths = {
        "#": 4, "PDB_ID": 18, "Score": 6, "Seq_id": 6, "Pubmed_id": 10,
        "Polymer": 10, "Assembly": 14, "Method": 18, "3D View": 16,
        "pH": 6, "Temp": 6, "Compounds": 57,
    }

    for (row, col), cell in tbl.get_celld().items():
        if row < 2:
            continue
        text_obj = cell.get_text()
        text = text_obj.get_text()
        col_name = sub_headers[col] if col < len(sub_headers) else None
        max_chars = wrap_widths.get(col_name, 20)
        if text:
            wrapped = "\n".join(textwrap.wrap(text, width=max_chars))
            text_obj.set_text(wrapped)
            text_obj.set_ha("center")
            text_obj.set_va("center")

    min_row_height = 0.06
    fig.canvas.draw()
    renderer = fig.canvas.get_renderer()
    row_heights = {}
    for (row, col), cell in tbl.get_celld().items():
        bbox = cell.get_text().get_window_extent(renderer=renderer)
        bbox_axes = bbox.transformed(ax.transAxes.inverted())
        height = max(bbox_axes.height * 1.7, min_row_height)
        row_heights[row] = max(row_heights.get(row, 0), height)

    for (row, col), cell in tbl.get_celld().items():
        cell.set_height(row_heights[row])
        cell.get_text().set_fontfamily("DejaVu Sans")

    tbl.auto_set_font_size(False)
    tbl.set_fontsize(11)
    tbl.scale(1, 1)

    cryst_cocktail_pdf = os.path.join(protein_dir, f"{protein_name}_Cryst_cocktail_Table.pdf")
    fig.savefig(cryst_cocktail_pdf, bbox_inches="tight")
    plt.close(fig)

    print(f"Table saved as: {cryst_cocktail_pdf}")
    print(f"PEG plot saved as: {peg_png}")
    print(f"Temp plot saved as: {temp_png}")

    return {
        "temp_png": temp_png,
        "peg_png": peg_png,
        "table_pdf": cryst_cocktail_pdf,
        "grouped_csv": grouped_csv_path,
    }
