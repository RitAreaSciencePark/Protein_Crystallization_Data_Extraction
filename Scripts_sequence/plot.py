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
        if not pd.isna(row["pH"]):
            return (row["pH"], 0.0, 0.0)
        if not pd.isna(row["pH_low"]) and not pd.isna(row["pH_high"]):
            ph = (row["pH_low"] + row["pH_high"]) / 2
            return (ph, ph - row["pH_low"], row["pH_high"] - ph)
        return (np.nan, 0.0, 0.0)

    df[["plot_pH", "err_low", "err_high"]] = df.apply(
        compute_plot_ph, axis=1, result_type="expand"
    )

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

    # Assign below-min pH / temp to "No pH" / "No temp" zones
    df["pH"] = df["plot_pH_numeric"].copy()
    df.loc[df["plot_pH_numeric"].isna() | (df["plot_pH_numeric"] < ph_min), "pH"] = no_ph_y

    df["temp"] = df["temp"]
    df.loc[df["temp"].isna() | (df["temp"] < temp_min), "temp"] = no_temp_x


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

    def group_by_condition(df, max_entries=10):

        condition_cols = ["temp", "pH", "method", "pubmed_id"]

        grouped = (
            df.groupby(condition_cols, dropna=False, as_index=False)
            .agg({
                "PDB_ID": lambda x: ", ".join(x),
                "pdbx_details": lambda x: " | ".join(x.fillna("")),
                "score": "mean",          # ✅ ADD THIS
            })
        )

        return grouped.head(max_entries)



    # ========================================================
    # 1️⃣ FULL PLOT
    # ========================================================
    fig, ax = plt.subplots(figsize=(14, 8))

    for _, row in df.iterrows():

        x = row["temp"]
        y = row["pH"]

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
    first10 = df.head(10).copy()
    first10["pdbx_details"] = first10["pdbx_details"].fillna("").astype(str)
    first10["marker"] = first10["method"].apply(assign_marker)
    
    first10_raw = df.head(50).copy()  # take more in case some merge reduces total
    first10_grouped = group_by_condition(first10_raw, max_entries=10)

    fig2 = plt.figure(figsize=(12, 10))
    gs = fig2.add_gridspec(2, 1, height_ratios=[1.5, 1.0], hspace=0.1)

    ax10 = fig2.add_subplot(gs[0])
    ax_table = fig2.add_subplot(gs[1])
    ax_table.axis("off")

    for _, row in first10_grouped.iterrows():
        x = row["temp"]
        y = row["pH"]
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
    table_cols = ["PDB_ID", "score", "pubmed_id", "temp", "pH", "pdbx_details", "method"]
    table_df = first10_grouped.reindex(columns=table_cols).fillna("")
    table_df["pdbx_details"] = table_df["pdbx_details"].apply(wrap_text) #wrap details
    table_df["score"] = table_df["score"].astype(float).apply(lambda x: f"{x:.3f}") #format score to 3 decimal

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
