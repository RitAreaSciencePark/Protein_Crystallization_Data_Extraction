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
    df["method"] = (
        df["method"]
        .fillna("unspecified")
        .astype(str)
        .str.lower()
        .str.replace(r"[^\w\s]", "", regex=True)
        .str.replace(r"\s+", " ", regex=True)
        .str.strip()
    )
    df.loc[df["method"] == "", "method"] = "unspecified"

    # --------------------------------------------------------
    # Numeric conversion
    # --------------------------------------------------------
    df["pH"] = pd.to_numeric(df["pH"], errors="coerce")
    df["temp"] = pd.to_numeric(df["temp"], errors="coerce")
    df["score"] = pd.to_numeric(df["score"], errors="coerce")
    df["score"] = df["score"].clip(0.5, 1)

    # --------------------------------------------------------
    # Parse pH range
    # --------------------------------------------------------
    def parse_ph_range(ph_range):
        if pd.isna(ph_range):
            return (np.nan, np.nan)
        m = re.search(r"([0-9.]+)\s*-\s*([0-9.]+)", str(ph_range))
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
    df["has_temp"] = ~df["temp"].isna()

    valid_ph = df["plot_pH"].dropna()
    valid_temp = df["temp"].dropna()

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
    df["temp"] = df["temp"]           # numeric temp from CSV or NaN

    # Plotting pH/temperature (with fallback for missing/too-low values)
    df["pH_plot"] = df["plot_pH_numeric"].copy()
    df.loc[df["plot_pH_numeric"].isna() | (df["plot_pH_numeric"] < ph_min), "pH_plot"] = no_ph_y

    df["temp_plot"] = df["temp"].copy()
    df.loc[df["temp"].isna() | (df["temp"] < temp_min), "temp_plot"] = no_temp_x

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

    def assign_marker(method):
        method = str(method).lower().strip()
        if method in method_marker_map:
            return method_marker_map[method]
        if method not in used_fallback:
            used_fallback[method] = len(used_fallback) % len(fallback_markers)
        index = used_fallback[method]
        return fallback_markers[index]

    # Apply to dataframe
    df["marker"] = df["method"].apply(assign_marker)

    cmap = plt.cm.viridis
    norm = plt.Normalize(0.5, 1)

    # -------------------------
    # Group top 10 conditions
    # -------------------------
    def first10_conditions_with_merged_pdb(df):
        """
        Return top 10 unique conditions based on pubmed_id, method, plot_pH_numeric, and COMPOUNDS (con_unit=mM).
        Merge PDB_IDs sharing the same condition while keeping all other columns intact.
        """
        cond_cols = ["pubmed_id", "method", "plot_pH_numeric", "COMPOUNDS (con_unit=mM)"]
        df_filled = df.copy()
        for col in cond_cols:
            df_filled[col] = df_filled[col].fillna("").astype(str)

        df_filled["_condition_key"] = df_filled[cond_cols].agg("||".join, axis=1)

        merged_rows = []
        for _, group in df_filled.groupby("_condition_key"):
            merged_pdb = ", ".join(group["PDB_ID"].astype(str).tolist())
            row = group.iloc[0].copy()  # take first row as representative
            row["PDB_ID"] = merged_pdb
            merged_rows.append(row)

        merged_df = pd.DataFrame(merged_rows)
        merged_df = merged_df.sort_values("score", ascending=False)
        top10 = merged_df.head(10).copy()
        top10.drop(columns=["_condition_key"], inplace=True)
        return top10

    # ========================================================
    # 1️⃣ FULL PLOT
    # ========================================================
    fig, ax = plt.subplots(figsize=(14, 8))

    for _, row in df.iterrows():
        x = row["temp_plot"] if "temp_plot" in row else row["temp"]
        y = row["pH_plot"] if "pH_plot" in row else row["pH"]
        ax.errorbar(x, y,
                    yerr=[[row["err_low"]], [row["err_high"]]] if row["has_ph"] else None,
                    fmt=row["marker"],
                    color=cmap(norm(row["score"])),
                    ecolor=cmap(norm(row["score"])),
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
    

    unique_methods = df["method"].unique()
    legend_items = []
    used_markers = {}
    for method in unique_methods:
        marker = assign_marker(method)
        if marker not in used_markers:
            legend_items.append(
                mlines.Line2D([], [], color="black", marker=marker, linestyle="None",
                            markersize=8, label=method.title())
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
    first10_grouped = first10_conditions_with_merged_pdb(df)  # your function

    fig, ax = plt.subplots(figsize=(14, 6))
    ax.axis("off")

    # 1️⃣ Table body rows
    table_rows = []
    for _, group in first10_grouped.iterrows():
        table_rows.append([
            group["PDB_ID"],
            f"{group['score']:.3f}" if not pd.isna(group["score"]) else "",
            group["pubmed_id"] if not pd.isna(group["pubmed_id"]) else "",
            group.get("Assembly", ""),
            group["plot_pH_numeric"] if not pd.isna(group["plot_pH_numeric"]) else "",
            group["temp"] if not pd.isna(group["temp"]) else "",
            group["method"],
            group["ligands"] if not pd.isna(group["ligands"]) else "",
            group["COMPOUNDS (con_unit=mM)"] if not pd.isna(group["COMPOUNDS (con_unit=mM)"]) else "",
        ])

    # 2️⃣ Headers
    main_header = ["", "", "", "", "", "", "CRYSTALLIZATION COCKTAILS", "", ""] 
    sub_headers = ["PDB_ID", "score", "pubmed_id", "Assembly","pH", "temp", "method", "ligands", "COMPOUNDS (con_unit=mM)"]

    # 3️⃣ Column widths
    col_widths = [0.09, 0.04, 0.07, 0.08, 0.04, 0.04, 0.10, 0.11, 0.35]

    # 4️⃣ Create table
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

    # 5️⃣ Experimental Conditions spanning pH → COMPOUNDS
    exp_start = sub_headers.index("pH")
    exp_end = sub_headers.index("COMPOUNDS (con_unit=mM)")

    # Keep only outer borders
    for c in range(exp_start, exp_end + 1):
        cell = tbl[0, c]
        if c == exp_start:
            cell.visible_edges = "LTB"  # left, top, bottom
        elif c == exp_end:
            cell.visible_edges = "RTB"  # right, top, bottom
        else:
            cell.visible_edges = "TB"   # top & bottom only
        cell.set_linewidth(1.2)

    # 7️⃣ Wrap data rows (row 2+)
    max_chars_dict = {
        "PDB_ID": 11,
        "method": 16,
        "ligands": 13,
        "COMPOUNDS (con_unit=mM)": 60}


    for r in range(len(table_rows)):
        for c, col_name in enumerate(sub_headers):

            cell = tbl[r + 2, c]
            text_obj = cell.get_text()

            # Wrap text
            text = text_obj.get_text()
            max_chars = max_chars_dict.get(col_name, 60)
            if len(text) > max_chars:
                wrapped_text = "\n".join(textwrap.wrap(text, width=max_chars))
                text_obj.set_text(wrapped_text)
                # Centered cells
                x, y = 0.5, 0.5
                text_obj.set_position((x, y))
                text_obj.set_ha("center")
                text_obj.set_va("center")
                cell._loc = None


    # 8️⃣ Scale & font
    tbl.auto_set_font_size(False)
    tbl.set_fontsize(8.5)
    tbl.scale(1, 3.5)

    # 9️⃣ Save PDF
    first10_pdf = os.path.join(
        os.path.dirname(output_csv_file),
        f"{protein_name}_FIRST10_TABLE.pdf")

    fig.savefig(first10_pdf, bbox_inches="tight")
    plt.close(fig)

    print(f"✔ Table saved as: {first10_pdf}")
