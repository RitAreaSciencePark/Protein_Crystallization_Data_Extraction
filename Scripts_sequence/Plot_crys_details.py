import pandas as pd
import matplotlib.pyplot as plt
import matplotlib.lines as mlines
import matplotlib.patheffects as pe
import numpy as np
import re
import os

def run_plot(csv_file):
    # ---------- Load data ----------
    print(f"▶ Loading CSV: {csv_file}")
    df = pd.read_csv(csv_file)

    # ---------- Representative PDB: first row ----------
    rep_pdb = df.iloc[0]["PDB_ID"] if not df.empty else "N/A"

    # ---------- Clean numeric ----------
    df["pH"] = pd.to_numeric(df["pH"], errors="coerce")
    df["temp"] = pd.to_numeric(df["temp"], errors="coerce")
    df["score"] = pd.to_numeric(df["score"], errors="coerce")
    df = df.dropna(subset=["temp", "PDB_ID", "score", "method"])

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

    # ---------- Compute plot pH ----------
    def compute_plot_ph(row):
        if not pd.isna(row["pH"]):
            ph = row["pH"]
            lo = ph - row["pH_low"] if not pd.isna(row["pH_low"]) else 0
            hi = row["pH_high"] - ph if not pd.isna(row["pH_high"]) else 0
        elif not pd.isna(row["pH_low"]) and not pd.isna(row["pH_high"]):
            ph = (row["pH_low"] + row["pH_high"]) / 2
            lo = ph - row["pH_low"]
            hi = row["pH_high"] - ph
        else:
            return pd.Series([np.nan, 0, 0])
        return pd.Series([ph, lo, hi])

    df[["plot_pH", "err_low", "err_high"]] = df.apply(compute_plot_ph, axis=1)
    df = df.dropna(subset=["plot_pH"])

    # ---------- Collapse identical conditions per method ----------
    group_cols = ["temp", "plot_pH", "method"]
    if "pubmed_id" in df.columns:
        group_cols.append("pubmed_id")

    grouped = (
        df.groupby(group_cols, as_index=False)
        .agg(
            {
                "PDB_ID": lambda x: sorted(set(x))[0],  # pick only one PDB per group
                "score": "mean",
                "err_low": "first",
                "err_high": "first",
            }
        )
    )

    # ---------- Marker selection ----------
    def choose_marker(method):
        if "hanging" in method.lower():
            return "^"
        elif "sitting" in method.lower():
            return "s"
        return "o"

    grouped["marker"] = grouped["method"].apply(choose_marker)

    # ---------- Colormap ----------
    cmap = plt.cm.viridis
    norm = plt.Normalize(grouped["score"].min(), grouped["score"].max())

    # ---------- Plot ----------
    fig, ax = plt.subplots(figsize=(10, 7))
    xlim = (grouped["temp"].min() - 5, grouped["temp"].max() + 5)
    ax.set_xlim(xlim)

    for _, row in grouped.iterrows():
        color = cmap(norm(row["score"]))

        # Plot point with error
        ax.errorbar(
            row["temp"],
            row["plot_pH"],
            yerr=[[row["err_low"]], [row["err_high"]]],
            fmt=row["marker"],
            color=color,
            ecolor=color,
            markersize=9,
            markeredgecolor="black",
            capsize=3,
            elinewidth=1,
        )

        # ---------- Label by default on the side ----------
        offset_x = 0.02 * (ax.get_xlim()[1] - ax.get_xlim()[0])
        ha = "left"
        if row["temp"] + offset_x > xlim[1]:
            offset_x = -0.02 * (ax.get_xlim()[1] - ax.get_xlim()[0])
            ha = "right"

        ax.text(
            row["temp"] + offset_x,
            row["plot_pH"],
            row["PDB_ID"],
            fontsize=8,
            color="black",
            ha=ha,
            va="center",
            path_effects=[pe.withStroke(linewidth=2, foreground="white")],
        )

    # ---------- Axes ----------
    ax.set_xlabel("Temperature (K)")
    ax.set_ylabel("pH")
    ax.set_title(f"pH vs Temperature (K) | PDB_ID = {rep_pdb}")
    ax.grid(False)

    # ---------- Method legend ----------
    legend_items = [
        mlines.Line2D([], [], marker="^", linestyle="None",
                      color="black", label="Hanging drop"),
        mlines.Line2D([], [], marker="s", linestyle="None",
                      color="black", label="Sitting drop"),
    ]

    ax.legend(
        handles=legend_items,
        title="Method",
        loc="lower center",
        ncol=2,
        bbox_to_anchor=(0.5, -0.25),
        frameon=False,
    )

    # ---------- Score colorbar ----------
    sm = plt.cm.ScalarMappable(cmap=cmap, norm=norm)
    sm.set_array([])
    cbar = fig.colorbar(sm, ax=ax, pad=0.02)
    cbar.set_label("Score", rotation=270, labelpad=15)

    # ---------- Save ----------
    out_pdf = os.path.join(
        os.path.dirname(csv_file),
        "pH_vs_Temperature_score_colorbar_labeled.pdf",
    )

    plt.tight_layout()
    plt.savefig(out_pdf)
    plt.show()

    print(f"Plot saved as: {out_pdf}")
