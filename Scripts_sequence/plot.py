import pandas as pd
import matplotlib.pyplot as plt
import matplotlib.lines as mlines
import matplotlib.patheffects as pe
import numpy as np
import re
import os
import textwrap


# ============================================================
# Helper: compute adaptive column widths
# ============================================================
def compute_col_widths(df, min_width=0.05, max_width=0.25, first_col_scale=1.7):
    widths = []
    for i, col in enumerate(df.columns):
        max_len = max(df[col].fillna("").astype(str).map(len).max(), len(col))
        w = min(max(min_width, max_len / 120), max_width)
        widths.append(w * first_col_scale if i == 0 else w)
    return widths


# ============================================================
# Helper: wrap long text
# ============================================================
def wrap_text(text, width=25, max_lines=2):
    if pd.isna(text) or text == "":
        return ""
    lines = textwrap.wrap(str(text), width=width)
    return "\n".join(lines[:max_lines])


# ============================================================
# MAIN FUNCTION
# ============================================================
def run_plot(csv_file):

    print(f"▶ Loading CSV: {csv_file}")
    df = pd.read_csv(csv_file)

    protein_name = os.path.basename(os.path.dirname(csv_file))

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
    def group_for_plotting(df):
        """
        Group by pubmed_id, temp, plot_pH_numeric, and method.
        All PDB_IDs with the same condition go into the same cell.
        Missing pH/temp are left empty in the table but assigned plotting positions.
        """
        condition_cols = ["pubmed_id", "temp", "plot_pH_numeric", "method"]

        def join_pdb_ids(series):
            clean = series.dropna()
            if clean.empty:
                return ""
            return ", ".join(dict.fromkeys(clean.astype(str)))  # preserve CSV order

        grouped = (
            df.groupby(condition_cols, dropna=False, as_index=False)
            .agg({
                "PDB_ID": join_pdb_ids,
                "pdbx_details": lambda x: " | ".join(x.dropna().astype(str)),
                "score": "mean"
            })
        )

        # Assign plotting positions for missing or below-min values
        grouped["temp_plot"] = grouped["temp"].copy()
        grouped.loc[grouped["temp_plot"].isna() | (grouped["temp_plot"] < temp_min), "temp_plot"] = no_temp_x

        grouped["pH_plot"] = grouped["plot_pH_numeric"].copy()
        grouped.loc[grouped["pH_plot"].isna() | (grouped["pH_plot"] < ph_min), "pH_plot"] = no_ph_y

        # Sort by mean score descending
        grouped = grouped.sort_values("score", ascending=False)

        return grouped.head(10)  # top 10 conditions


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
        os.path.dirname(csv_file),
        f"{protein_name}_FULL_PLOT.png",
    )
    fig.savefig(full_png, dpi=300, bbox_inches="tight")
    plt.close(fig)

    print(f"✔ Full plot saved as: {full_png}")

    # ========================================================
    # 2️⃣ FIRST 10 + TABLE
    # ========================================================
    
    # 1️⃣ Group for plotting
    first10_grouped = group_for_plotting(df)

    fig2 = plt.figure(figsize=(12, 10))
    gs = fig2.add_gridspec(2, 1, height_ratios=[1.5, 1.0], hspace=0.1)

    ax10 = fig2.add_subplot(gs[0])
    ax_table = fig2.add_subplot(gs[1])
    ax_table.axis("off")

    for _, row in first10_grouped.iterrows():
        x = row["temp_plot"] if "temp_plot" in row else row["temp"]
        y = row["pH_plot"] if "pH_plot" in row else row["pH"]
        marker = assign_marker(row["method"])
        c = cmap(norm(row["score"]))

        ax10.scatter( x, y, color=c, edgecolor="black", s=80, marker=marker)

        # Add PDB_ID label
        ax10.text(x, y, row["PDB_ID"], fontsize=8, path_effects=[pe.withStroke(linewidth=2, foreground="white")])


        ax10.set_title(
            f"First 10 PDB entries with different conditions\nProtein: {protein_name}",
            fontsize=12,
            fontweight="bold",
        )
    # Dashed lines for No pH / No Temp
    ax10.axhline(no_ph_y + 0.25, linestyle="--", color="black")
    ax10.axvline(no_temp_x + 2, linestyle="--", color="black")

    # Get existing ticks, Remove any tick below real minimum, Add artificial no_temp position, Create labels
    ax10.set_xlim(no_temp_x - 2, valid_temp.max() + 2)
    xticks = [no_temp_x] + list(range(temp_min_tick, temp_max_tick + 5, 10))
    ax10.set_xticks(xticks)
    ax10.set_xticklabels(["No Temp"] + [str(t) for t in xticks[1:]])


    # Get existing ticks, Remove any tick below real minimum, Add artificial no_pH position, Create labels
    ax10.set_ylim(no_ph_y - 0.25, valid_ph.max() + 0.25)
    yticks = [no_ph_y] + list(np.arange(ph_min_tick, ph_max_tick + 0.5, 1))
    ax10.set_yticks(yticks)
    ax10.set_yticklabels(["No pH"] + [f"{t:.1f}" for t in yticks[1:]])

    ax10.set_xlabel("Temperature (K)")
    ax10.set_ylabel("pH")

    # ---------- Colorbar ----------
    sm2 = plt.cm.ScalarMappable(cmap=cmap, norm=norm)
    sm2.set_array([])
    cbar2 = fig2.colorbar(sm2, ax=ax10, pad=0.02)
    cbar2.set_label("Score", rotation=270, labelpad=15)
    cbar2.set_ticks([0.5, 0.6, 0.7, 0.8, 0.9, 1.0])
    cbar2.set_ticklabels([f"{t:.1f}" for t in [0.5, 0.6, 0.7, 0.8, 0.9, 1.0]])

    # ---------- Table ----------
    def build_table(top_groups):
        """
        Expand PDB_IDs into individual rows for the table.
        Keep missing pH/temp as empty strings.
        """
        table_rows = []
        for _, group in top_groups.iterrows():
                table_rows.append({
                    "PDB_ID": group["PDB_ID"],  # <- use group["PDB_ID"]
                    "score": f"{group['score']:.3f}" if not pd.isna(group["score"]) else "",
                    "pubmed_id": group["pubmed_id"] if not pd.isna(group["pubmed_id"]) else "",
                    "temp": group["temp"] if not pd.isna(group["temp"]) else "",
                    "pH": group["plot_pH_numeric"] if not pd.isna(group["plot_pH_numeric"]) else "",
                    "pdbx_details": wrap_text(group["pdbx_details"]) if not pd.isna(group["pdbx_details"]) else "",
                    "method": group["method"]
                })

        table_df = pd.DataFrame(table_rows)
        # Sort table by score descending
        table_df = table_df.sort_values("score", ascending=False)
        return table_df


    # Build the table
    table_df = build_table(first10_grouped)
    col_widths = compute_col_widths(table_df)

    tbl = ax_table.table(
        cellText=table_df.values,
        colLabels=table_df.columns,
        colWidths=col_widths,
        cellLoc="left",
        colLoc="left",
        loc="upper center",
    )

    tbl.auto_set_font_size(False)
    tbl.set_fontsize(8.5)
    tbl.scale(1, 2.1)

    first10_png = os.path.join(
        os.path.dirname(csv_file),
        f"{protein_name}_FIRST10_TABLE.png",
    )
    fig2.savefig(first10_png, dpi=300, bbox_inches="tight")
    plt.close(fig2)

    print(f"✔ First 10 plot saved as: {first10_png}")
