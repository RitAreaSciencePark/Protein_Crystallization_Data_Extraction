import pandas as pd
import matplotlib.pyplot as plt
import matplotlib.lines as mlines
import numpy as np
import re
import os
def run_plot(csv_file):
    # ---------- Load data ----------
    print(f"▶ Loading CSV: {csv_file}")
    df = pd.read_csv(csv_file)

    # ---------- Clean numeric ----------
    df["pH"] = pd.to_numeric(df["pH"], errors="coerce")
    df["temp"] = pd.to_numeric(df["temp"], errors="coerce")
    df["score"] = pd.to_numeric(df["score"], errors="coerce")

    df = df.dropna(subset=["pH", "temp", "method", "PDB_ID", "score"])

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

    # ---------- Colormap for score ----------
    cmap = plt.cm.viridis
    norm = plt.Normalize(df["score"].min(), df["score"].max())

    # ---------- Plot ----------
    #plt.figure(figsize=(10, 7))
    fig, ax = plt.subplots(figsize=(10, 7))

    for _, row in df.iterrows():
        color = cmap(norm(row["score"]))

        ax.errorbar(
            row["temp"],
            row["pH"],
            fmt=row["marker"],
            color=color,
            ecolor=color,
            markersize=8,
            markeredgecolor="black",
            capsize=3,
            elinewidth=1
        )

        # ---------- Method legend only ----------
        method_legend = [
            mlines.Line2D([], [], color='black', marker='^',
                        linestyle='None', markersize=8, label='Hanging drop'),
            mlines.Line2D([], [], color='black', marker='s',
                        linestyle='None', markersize=8, label='Sitting drop')
        ]

        ax.legend(handles=method_legend, title="Method", loc="upper left")

        # ---------- PDB_ID at each point ----------
        ax.annotate(
            row["PDB_ID"],
            (row["temp"], row["pH"]),
            xytext=(5, 5),
            textcoords="offset points",
            fontsize=7,
            ha="left",
            va="bottom"
        )

    # ---------- Axes ----------
    ax.set_xlabel("Temperature (K)")
    ax.set_ylabel("pH")
    ax.set_title("pH vs Temperature (K)\nprotein = PLpro | PDB_ID = 9BRX")
    ax.grid(True)

    # ---------- X-axis ticks ----------
    unique_temps = sorted(df["temp"].unique())
    ax.set_xticks(unique_temps)
    ax.set_xticklabels([f"{t:.1f}" for t in unique_temps], rotation=45)

    # ---------- Method legend only ----------
    method_legend = [
        mlines.Line2D([], [], color='black', marker='^',
                    linestyle='None', markersize=8, label='Hanging drop'),
        mlines.Line2D([], [], color='black', marker='s',
                    linestyle='None', markersize=8, label='Sitting drop')
    ]

    ax.legend(
        handles=method_legend,
        title="Method",
        loc="lower center",
        ncol=2,
        bbox_to_anchor=(0.5, -0.23),
        frameon=False
    )
    # ---------- Score colorbar ----------
    sm = plt.cm.ScalarMappable(cmap=cmap, norm=norm)
    sm.set_array([])

    cbar = fig.colorbar(sm, ax=ax, pad=0.02)
    cbar.set_label("Score", rotation=270, labelpad=15)

    # ---------- Save ----------
    output_pdf = os.path.join(os.path.dirname(csv_file), "pH_vs_Temperature_score_colorbar_labeled.pdf")
    plt.tight_layout()
    plt.savefig(output_pdf)
    plt.show()

    print(f"Plot saved as: {output_pdf}")
