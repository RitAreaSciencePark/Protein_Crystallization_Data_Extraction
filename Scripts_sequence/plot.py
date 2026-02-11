import pandas as pd
import matplotlib.pyplot as plt
import matplotlib.lines as mlines
import matplotlib.patheffects as pe
import numpy as np
import re
import os
import textwrap
from matplotlib.backends.backend_pdf import PdfPages

# Helper function to calculate automatic column widths
def compute_col_widths(df, min_width=0.05, max_width=0.6, first_col_scale=1.7):
    widths = []
    for i, col in enumerate(df.columns):
        max_len = max(df[col].fillna("").astype(str).map(len).max(), len(col))
        w = min(max(min_width, max_len / 120), max_width)
        widths.append(w * first_col_scale if i == 0 else w)
    return widths


# Wrap text to max 2 lines for table
def wrap_text(text, width=40, max_lines=2):
    if pd.isna(text) or text == "":
        return ""
    lines = textwrap.wrap(str(text), width=width)
    return "\n".join(lines[:max_lines])

def run_plot(csv_file):
    print(f"▶ Loading CSV: {csv_file}")
    df = pd.read_csv(csv_file)

    df["method"] = (
    df["method"]
    .fillna("unspecified")
    .astype(str)
    .str.lower()
    .str.replace(r"[^\w\s]", "", regex=True)  # remove ?, ', ., commas, etc.
    .str.replace(r"\s+", " ", regex=True)     # normalize spaces
    .str.strip())

    #............replace empty strings after cleaning
    df.loc[df["method"] == "", "method"] = "unspecified"


    # ---------- Representative PDB ----------
    rep_pdb = df.iloc[0]["PDB_ID"] if not df.empty else "N/A"

    # ---------- Table: entries missing pH or temperature ----------
    missing_mask = df["pH"].isna() | df["temp"].isna()
    missing_table = df.loc[
        missing_mask,
        ["PDB_ID", "score", "pubmed_id", "pH", "method", "temp", "pdbx_details"]
    ].copy()
    missing_table["pdbx_details"] = missing_table["pdbx_details"].fillna("").astype(str)

    print(f"Missing table rows: {missing_table.shape[0]}")

    # ---------- Clean numeric ----------
    df["pH"] = pd.to_numeric(df["pH"], errors="coerce")
    df["temp"] = pd.to_numeric(df["temp"], errors="coerce")
    df["score"] = pd.to_numeric(df["score"], errors="coerce")
    df = df.dropna(subset=["temp", "PDB_ID", "score"])

    # ---------- Keep hanging / sitting / missing ----------
    df["method"] = df["method"].fillna("unspecified").astype(str).str.strip().str.lower()
    df["method"] = df["method"].str.replace("vapor diffusion", "", regex=False).str.strip()
    unique_methods = sorted(df["method"].unique())
    print("Detected methods:", unique_methods)


    # ---------- Parse pH range ----------
    def parse_ph_range(ph_range):
        if pd.isna(ph_range):
            return None, None
        m = re.match(r"([0-9.]+)\s*-\s*([0-9.]+)", str(ph_range))
        if m:
            return float(m.group(1)), float(m.group(2))
        return None, None

    df[["pH_low", "pH_high"]] = df["pdbx_pH_range"].apply(lambda x: pd.Series(parse_ph_range(x)))

    # ---------- Compute plot pH ----------
    def compute_plot_ph(row):
        if not pd.isna(row["pH"]):
            return float(row["pH"]), 0, 0
        if not pd.isna(row["pH_low"]) and not pd.isna(row["pH_high"]):
            ph = (row["pH_low"] + row["pH_high"]) / 2
            return float(ph), float(ph - row["pH_low"]), float(row["pH_high"] - ph)
        return np.nan, 0, 0

    df[["plot_pH", "err_low", "err_high"]] = df.apply(lambda r: pd.Series(compute_plot_ph(r)), axis=1)
    df = df.dropna(subset=["plot_pH"])

    # ---------- Marker selection ----------
    # ---------- Marker selection ----------
    marker_cycle = ["o", "^", "s", "D", "v", "P", "X", "*", "<", ">"]

    method_to_marker = {method: marker_cycle[i % len(marker_cycle)]
        for i, method in enumerate(unique_methods)}

    # enforce unspecified as circle
    method_to_marker["unspecified"] = "o"

    df["marker"] = df["method"].map(method_to_marker)

    # ---------- Color map ----------
    cmap = plt.cm.viridis
    norm = plt.Normalize(df["score"].min(), df["score"].max())

    # ---------- First 10 PDBs ----------
    first10 = df.head(10).copy()
    first10["marker"] = first10["method"].map(method_to_marker)
    first10["pdbx_details"] = first10["pdbx_details"].fillna("").astype(str)

    # ---------- Layout ----------
    fig = plt.figure(figsize=(22, 10))
    gs = fig.add_gridspec(1, 2, width_ratios=[2, 1.0], wspace=0.1)
    ax_main = fig.add_subplot(gs[0, 0])

    right = gs[0, 1].subgridspec(2, 1, height_ratios=[1.5, 1.0], hspace=0.1)
    ax10 = fig.add_subplot(right[0])
    ax10_tbl = fig.add_subplot(right[1])
    ax10_tbl.axis("off")

    # ---------- Main plot ----------
    for _, row in df.iterrows():
        c = cmap(norm(row["score"]))
        ax_main.errorbar(
            row["temp"],
            row["plot_pH"],
            yerr=[[row["err_low"]], [row["err_high"]]],
            fmt=row["marker"],
            color=c,
            ecolor=c,
            markeredgecolor="black",
            markersize=9,
            capsize=3,
        )
        ax_main.text(
            row["temp"],
            row["plot_pH"],
            row["PDB_ID"],
            fontsize=8,
            color="black",
            path_effects=[pe.withStroke(linewidth=2, foreground="white")],
        )

    ax_main.set_xlabel("Temperature (K)")
    ax_main.set_ylabel("pH")
    ax_main.set_title(f"pH vs Temperature (K) | PDB_ID = {rep_pdb}")

    # ---------- Legend ----------
    legend_items =[]
    for method, marker in method_to_marker.items():
        label = method.title()
        legend_items.append(
        mlines.Line2D(
            [],
            [],
            marker=marker,
            linestyle="None",
            markeredgecolor="black",
            color="black",
            label=label,
            markersize=9,
        )
    )
    ax_main.legend(handles=legend_items, title="Method", loc="lower center", ncol=3, bbox_to_anchor=(0.5, -0.3), frameon=False)

    # ---------- Score colorbar ----------
    sm = plt.cm.ScalarMappable(cmap=cmap, norm=norm)
    sm.set_array([])
    cbar = fig.colorbar(sm, ax=ax_main, pad=0.02)
    cbar.set_label("Score", rotation=270, labelpad=15)

    # ---------- Right figure: first 10 PDBs with all labels ----------
    y_offsets = {}
    for _, row in first10.iterrows():
        c = cmap(norm(row["score"]))
        x = row["temp"]
        y = row["plot_pH"]
        ax10.errorbar(
            x, y,
            fmt=row["marker"],
            color=c,
            markeredgecolor="black",
            markersize=8,
        )

        # Handle overlapping labels
        y_key = round(y, 3)
        if y_key not in y_offsets:
            y_offsets[y_key] = 0.05
        else:
            y_offsets[y_key] += 0.15  # shift up for overlap

        y_label = y + y_offsets[y_key]
        # prevent leaving plot area
        y_label = max(min(y_label, ax10.get_ylim()[1]), ax10.get_ylim()[0])
        ax10.text(
            x, y_label,
            row["PDB_ID"],
            fontsize=8,
            ha="center",
        )

    ax10.set_title("First 10 PDB entries (individual points)")
    ax10.set_xlabel("Temperature (K)")
    ax10.set_ylabel("pH")

    # ---------- Table below right figure ----------
    table_df = first10[["PDB_ID", "score", "temp", "pH", "pdbx_pH_range", "pdbx_details"]].fillna("")
    table_df["pdbx_details"] = table_df["pdbx_details"].apply(lambda x: wrap_text(x, width=40, max_lines=2))
    col_widths_10 = compute_col_widths(table_df, first_col_scale=1.7)

    tbl = ax10_tbl.table(
        cellText=table_df.values,
        colLabels=table_df.columns,
        cellLoc="left",
        colLoc="left",
        colWidths=col_widths_10,
        loc="upper center",
    )
    tbl.auto_set_font_size(False)
    tbl.set_fontsize(8.5)
    tbl.scale(1, 2.2)
    for (r, c), cell in tbl.get_celld().items():
        cell.set_edgecolor("black")
        if r == 0:
            cell.set_text_props(weight="bold")
            cell.set_facecolor("#f0f0f0")

    # ---------- Save PDF ----------
    out_pdf = os.path.join(os.path.dirname(csv_file), "pH_vs_Temperature_first10.pdf")
    with PdfPages(out_pdf) as pdf:
        pdf.savefig(fig, bbox_inches="tight")

        # ---------- Page 2: table of entries missing pH or temperature ----------
        if not missing_table.empty:
            fig2, ax2 = plt.subplots(figsize=(16, 11))
            ax2.axis("off")
            ax2.set_title("PDB entries missing pH or Temperature", fontsize=16, pad=25)
            col_widths_missing = compute_col_widths(missing_table)
            t = ax2.table(
                cellText=missing_table.fillna("").values,
                colLabels=missing_table.columns,
                cellLoc="left",
                colLoc="left",
                colWidths=col_widths_missing,
                loc="upper center",
            )
            t.auto_set_font_size(False)
            t.set_fontsize(9)
            t.scale(1, 2.0)
            for (r, c), cell in t.get_celld().items():
                cell.set_edgecolor("black")
                if r == 0:
                    cell.set_text_props(weight="bold")
                    cell.set_facecolor("#f0f0f0")
            pdf.savefig(fig2, bbox_inches="tight")
            plt.close(fig2)

    plt.show()
    print(f"Plot saved as: {out_pdf}")

