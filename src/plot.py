import pandas as pd
import matplotlib.pyplot as plt
import matplotlib.lines as mlines
import matplotlib.patheffects as pe
import numpy as np
import re
import os
import textwrap

# Automatically compute approximate column widths based on max text length
def compute_col_widths(df, scale=0.01):
    """
    Compute column widths for matplotlib table based on text length.
    scale: adjust to make table wider/narrower
    """
    widths = {}
    for col in df.columns:
        max_len = max(df[col].astype(str).apply(len).max(), len(col))
        widths[col] = max_len * scale  # simple scaling factor
    return widths

# ============================================================
# MAIN FUNCTION
# ============================================================
def run_plot(output_csv_file):

    print(f"▶ Loading CSV: {output_csv_file}")
    df = pd.read_csv(output_csv_file)
    protein_name = os.path.basename(output_csv_file).split("_crystallization_data")[0]

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
        
        # Match numbers with optional decimals
        m = re.search(r"(\d+(?:\.\d+)?)\s*[-–]\s*(\d+(?:\.\d+)?)", str(ph_range))
        if m:
            return (float(m.group(1)), float(m.group(2)))
        
        return (np.nan, np.nan)

    if "pdbx_pH_range" in df.columns:
        ranges = df["pdbx_pH_range"].apply(parse_ph_range)
        df["pH_low"] = ranges.apply(lambda x: x[0])
        df["pH_high"] = ranges.apply(lambda x: x[1])
    else:
        df["pH_low"] = np.nan
        df["pH_high"] = np.nan

    # --------------------------------------------------------
    # Compute plotting pH
    # --------------------------------------------------------
    def compute_plot_ph(row):
    # Case 1: exact pH value present
        if pd.notna(row.get("pH")):
            return pd.Series([row["pH"], 0.0, 0.0])

        # Case 2: pH range available
        if pd.notna(row.get("pH_low")) and pd.notna(row.get("pH_high")):
            ph_low = row["pH_low"]
            ph_high = row["pH_high"]
            ph_mid = (ph_low + ph_high) / 2
            return pd.Series([ph_mid, ph_mid - ph_low, ph_high - ph_mid])
        # Case 3: no pH info
        return pd.Series([np.nan, 0.0, 0.0])
        
    df[["plot_pH", "err_low", "err_high"]] = df.apply(compute_plot_ph, axis=1, result_type="expand")

    df["has_ph"] = ~df["plot_pH"].isna()
    df["has_temp"] = ~df["Temp"].isna()

    valid_ph = df["plot_pH"].dropna()
    valid_temp = df["Temp"].dropna()

    ph_min = 4.5
    ph_max = valid_ph.max()
    temp_min = 275
    temp_max = valid_temp.max()

    # Define compact "No Temp" /"No pH" region
    no_ph_y = ph_min - 0.5
    no_temp_x = temp_min - 4
    

    # Compute numeric pH / temp for plotting
    df["plot_pH_numeric"] = df.apply(lambda row: row["plot_pH"][0] if isinstance(row["plot_pH"], tuple) else row["plot_pH"], axis=1)

    # Keep original pH and temp for table/reporting
    df["pH"] = df["plot_pH_numeric"]  # numeric pH from CSV or NaN
    df["Temp"] = df["Temp"]           # numeric temp from CSV or NaN

    # Plotting pH/temperature (with fallback for missing/too-low values)
    df["pH_plot"] = df["plot_pH_numeric"].copy()
    df.loc[df["plot_pH_numeric"].isna() | (df["plot_pH_numeric"] < ph_min), "pH_plot"] = no_ph_y

    df["temp_plot"] = df["Temp"].copy()
    df.loc[df["Temp"].isna() | (df["Temp"] < temp_min), "temp_plot"] = no_temp_x

    # Round real temperature scale to 5K grid
    temp_min_tick = int(np.floor(temp_min / 5) * 5)
    temp_max_tick = int(np.ceil(temp_max / 5) * 5)

    # Round bounds nicely
    ph_min_tick = np.floor(ph_min * 2) / 2     # round down to nearest 0.5
    ph_max_tick = np.ceil(ph_max * 2) / 2      # round up to nearest 0.5


    # --------------------------------------------------------
    # Marker assignment
    # --------------------------------------------------------
        # Standardized method → marker mapping
    method_marker_map = {
        "Vapor Diffusion Hanging Drop": "s",  # square
        "Hanging Drop": "s",
        "Hanging Drop Vapor Diffusion": "s",
        "Vapor Diffusion Sitting Drop": "^",  # triangle
        "Sitting Drop": "^",
        "Sitting Drop Vapor Diffusion": "^",
        "vapor diffusion": "o",               # circle
        "unspecified": "X",
        "": "X",
    }

    # Fallback markers for other methods
    fallback_markers = ["D", "v", "P", "*", "<", ">"]
    used_fallback = {}  # track used methods

    def assign_marker(Method):
        Method = str(Method).lower().strip()
        if Method in method_marker_map:
            return method_marker_map[Method]
        if Method not in used_fallback:
            used_fallback[Method] = len(used_fallback) % len(fallback_markers)
        index = used_fallback[Method]
        return fallback_markers[index]

    # Apply to dataframe
    df["marker"] = df["Method"].apply(assign_marker)

    cmap = plt.cm.viridis
    norm = plt.Normalize(0.5, 1)

    # -------------------------
    # Group top 10 conditions
    # 
    def all_conditions_with_merged_pdb(df):
        """
        Return ALL unique conditions based on pubmed_id, method, plot_pH_numeric, 
        and COMPOUNDS (con_unit=mM).
        Merge PDB_IDs sharing the same condition while keeping all other columns intact.
        """
        cond_cols = ["Pubmed_id", "Method", "plot_pH_numeric", "Compounds(con_unit=mM)"]
        
        df_filled = df.copy()
        for col in cond_cols:
            df_filled[col] = df_filled[col].fillna("").astype(str)

        df_filled["_condition_key"] = df_filled[cond_cols].agg("||".join, axis=1)

        merged_rows = []
        for _, group in df_filled.groupby("_condition_key"):
            merged_pdb = ", ".join(group["PDB_ID"].astype(str).tolist())
            row = group.iloc[0].copy()  # representative row
            row["PDB_ID"] = merged_pdb
            merged_rows.append(row)

        merged_df = pd.DataFrame(merged_rows)
        merged_df = merged_df.sort_values("Score", ascending=False)

        # return everything instead of top 10
        merged_df.drop(columns=["_condition_key"], inplace=True)
        return merged_df
    # ========================================================
    # 1️⃣ FULL PLOT
    # ========================================================
    fig, ax = plt.subplots(figsize=(14, 8))

    for _, row in df.iterrows():
        x = row["temp_plot"] if "temp_plot" in row else row["Temp"]
        y = row["pH_plot"] if "pH_plot" in row else row["pH"]
        ax.errorbar(x, y,
                    yerr=[[row["err_low"]], [row["err_high"]]] if row["has_ph"] else None,
                    fmt=row["marker"],
                    color=cmap(norm(row["Score"])),
                    ecolor=cmap(norm(row["Score"])),
                    markeredgecolor="black",
                    markersize=9,
                    capsize=3)

        ax.text( x, y, row["PDB_ID"], fontsize=8, path_effects=[pe.withStroke(linewidth=2, foreground="white")],)

    ax.axhline(no_ph_y + 0.25, linestyle="--", color="black")
    ax.axvline(no_temp_x + 2, linestyle="--", color="black")

    # Get existing ticks, Remove any tick below real minimum, Add artificial no_temp position, Create labels Build consistent 2K temperature scale
    ax.set_xlim(no_temp_x - 2, valid_temp.max() + 5)
    xticks = [no_temp_x] + list(range(temp_min_tick, temp_max_tick + 5, 10))
    ax.set_xticks(xticks)
    ax.set_xticklabels(["No Temp"] + [str(t) for t in xticks[1:]])
   
    # Get existing ticks, Remove any tick below real minimum, Add artificial no_pH position, Create labels
    ax.set_ylim(no_ph_y - 0.25, valid_ph.max() + 0.25)
    yticks = [no_ph_y] + list(np.arange(ph_min_tick, ph_max_tick + 0.5, 1))
    ax.set_yticks(yticks)
    ax.set_yticklabels(["No pH"] + [f"{t:.1f}" for t in yticks[1:]])

    ax.set_xlabel("Temperature (K)")
    ax.set_ylabel("pH")
    ax.set_title(
        f"pH vs Temperature (K)\nProtein: {protein_name}",
        fontsize=14,
        fontweight="bold",
    )
    

    unique_methods = df["Method"].unique()
    legend_items = []
    used_markers = {}
    for Method in unique_methods:
        marker = assign_marker(Method)
        if marker not in used_markers:
            legend_items.append(
                mlines.Line2D([], [], color="black", marker=marker, linestyle="None",
                            markersize=8, label=Method.title())
            )
            used_markers[marker] = True

    # Place legend **just above the x-axis label**
    ax.legend(handles=legend_items,
            title="Method",
            loc='lower center',        # relative to axes
            bbox_to_anchor=(0.5, -0.25),  # slightly below x-axis (adjust as needed)
            ncol=3,
            frameon=False)

    # Fixed colorbar
    sm = plt.cm.ScalarMappable(cmap=cmap, norm=norm)
    sm.set_array([])
    cbar = fig.colorbar(sm, ax=ax, pad=0.02)
    cbar.set_label("Score", rotation=270, labelpad=15)
    fixed_ticks = [0.5, 0.6, 0.7, 0.8, 0.9, 1.0]
    cbar.set_ticks(fixed_ticks)
    cbar.set_ticklabels([f"{t:.1f}" for t in fixed_ticks])

    full_png = os.path.join(
        os.path.dirname(output_csv_file),
        f"{protein_name}_FULL_PLOT.png",)
    fig.savefig(full_png, dpi=300, bbox_inches="tight")
    plt.close(fig)
    print(f"✔ Full plot saved as: {full_png}")

    # ========================================================
    # 2️⃣ FIRST 10 + TABLE
    # ========================================================
    # Example: first10_grouped DataFrame
    first10_grouped = all_conditions_with_merged_pdb(df)  # your function

    fig, ax = plt.subplots(figsize=(16, 10))
    ax.axis("off")

    # Table body rows
    table_rows = []
    for _, group in first10_grouped.iterrows():
        table_rows.append([
            group["PDB_ID"],
            f"{group['Score']:.3f}" if not pd.isna(group["Score"]) else "",  
            f"{group['Seq_id']:.3f}" if not pd.isna(group["Seq_id"]) else "", 
            group["Pubmed_id"] if not pd.isna(group["Pubmed_id"]) else "",  
            group.get("Polymer", ""), group.get("Assembly", ""), group["Method"], group["Ligands"] if not pd.isna(group["Ligands"]) else "",
            group["plot_pH_numeric"] if not pd.isna(group["plot_pH_numeric"]) else "", group["Temp"] if not pd.isna(group["Temp"]) else "",
            group["Compounds(con_unit=mM)"] if not pd.isna(group["Compounds(con_unit=mM)"]) else "",])

    # 2️⃣ Headers
    main_header = ["", "", "", "", "", "", "", "", "", "", "CRYSTALLIZATION COCKTAILS"] 
    sub_headers = ["PDB_ID", "Score", "Seq_id", "Pubmed_id", "Polymer", "Assembly","Method", "Ligands","pH", "Temp", "Compounds(con_unit=mM)"]

    # 3️⃣ Column widths
    col_widths = [0.16, 0.05, 0.05, 0.08, 0.07, 0.09, 0.11, 0.10, 0.04, 0.04, 0.36]

    # Create table
    tbl = ax.table(
        cellText=[main_header, sub_headers] + table_rows,
        colWidths=col_widths,
        cellLoc="center",
        loc="center")

    # Main header font (row 0)
    for c in range(len(main_header)):
        cell = tbl[0, c]
        txt = cell.get_text()
        txt.set_fontsize(18)    
        txt.set_weight("bold")
        

    for c in range(len(sub_headers)):
        cell = tbl[1, c]
        txt = cell.get_text()
        txt.set_weight("bold")
        txt.set_fontsize(12)   

    # 🎨 Define a color per column
    col_colors = {
        "PDB_ID": "#98ACBA",
        "Score": "#B7EFBC",
        "Pubmed_id": "#D8BC8E",
        "Assembly": "#CA94D2",
        "Method": "#9AD7DF",
        "Ligands": "#E9A6BC",
        "pH": "#818FD8",
        "Temp": "#CBC8B1",
        "Compounds(con_unit=mM)": "#A7F054",
    }

    # Apply colors to columns (from sub_header row down)
    for col, col_name in enumerate(sub_headers):
        color = col_colors.get(col_name, "white")

        for row in range(1, len(table_rows) + 2):  # row 1 = sub_header, rest = data
            cell = tbl[row, col]
            cell.set_facecolor(color)

    # Experimental Conditions spanning pH → COMPOUNDS
    exp_start = sub_headers.index("pH")
    exp_end = sub_headers.index("Compounds(con_unit=mM)")

    # Keep only outer borders
    group_color = "#F2EEED"  # soft green (change as you like)

    for c in range(exp_start, exp_end + 1):
        cell = tbl[0, c]

        # 🔹 Apply DIFFERENT color for grouped header
        cell.set_facecolor(group_color)

        # 🔹 Keep your border logic
        if c == exp_start:
            cell.visible_edges = "LTB"
        elif c == exp_end:
            cell.visible_edges = "RTB"
        else:
            cell.visible_edges = "TB"

        cell.set_linewidth(1.2)

    # Remove all vertical separators from PDB_IDs → ligands

    pdb_start = sub_headers.index("PDB_ID")
    pdb_end = sub_headers.index("Ligands")

    for c in range(pdb_start, pdb_end + 1):
        cell = tbl[0, c]
        if c == pdb_end:
            cell.visible_edges = "RB"  # right, bottom
        else:
            cell.visible_edges = "B"   #  bottom only (no verticals)
        cell.set_linewidth(1.2)

    # Wrap text consistently per column (stable layout)
    wrap_widths = {
        "PDB_ID": 18,
        "Score": 6,
        "Seq_id": 6,
        "Pubmed_id": 10,
        "Polymer": 10,
        "Assembly": 14,
        "Method": 18,
        "Ligands": 15,
        "pH": 6,
        "Temp": 6,
        "Compounds(con_unit=mM)": 57,}

    for (row, col), cell in tbl.get_celld().items():

        if row < 2:   #  skip main_header (0) and sub_headers (1)
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

    # 🔹 Minimum row height in axes fraction
    min_row_height = 0.06  # adjust as needed (0.03 = ~3% of axes height)
    fig.canvas.draw()
    renderer = fig.canvas.get_renderer()

    row_heights = {}

    for (row, col), cell in tbl.get_celld().items():
        text_obj = cell.get_text()
        bbox = text_obj.get_window_extent(renderer=renderer)

        bbox_axes = bbox.transformed(ax.transAxes.inverted())
        height = bbox_axes.height * 1.7 # padding

        height = max(height, min_row_height)

        row_heights[row] = max(row_heights.get(row, 0), height)

           # Scale & font
    for (row, col), cell in tbl.get_celld().items():
        cell.set_height(row_heights[row])
        cell.get_text().set_fontfamily("DejaVu Sans")
    
    tbl.auto_set_font_size(False)
    tbl.set_fontsize(11)
    tbl.scale(1, 1)

    # Save PDF
    Cryst_cocktail = os.path.join(
        os.path.dirname(output_csv_file),
        f"{protein_name}_Cryst_cocktail_Table.pdf")

    fig.savefig(Cryst_cocktail, bbox_inches="tight")
    plt.close(fig)

    print(f"✔ Table saved as: {Cryst_cocktail}")
