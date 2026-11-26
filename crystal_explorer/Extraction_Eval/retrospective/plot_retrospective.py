"""Figure for the retrospective section, from retrospective_summary.json and
context_summary.json (context_analysis.py).

    python plot_retrospective.py   # -> ../../../paper/figures/retrospective.pdf (+ .png)
"""
import json
import os

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt  # noqa: E402
import numpy as np  # noqa: E402

HERE = os.path.dirname(os.path.abspath(__file__))
OUT = os.path.join(HERE, "..", "..", "..", "paper", "figures", "retrospective")

# reference palette, first two categorical slots (validated all-pairs)
BLUE, ORANGE, GREEN = "#2a78d6", "#eb6834", "#1f9e6e"
INK, INK2, GRID = "#0b0b0b", "#52514e", "#e4e3df"

plt.rcParams.update({
    "font.family": "serif", "font.size": 7.5, "axes.edgecolor": INK2,
    "axes.labelcolor": INK, "xtick.color": INK2, "ytick.color": INK2,
    "axes.spines.top": False, "axes.spines.right": False, "axes.linewidth": 0.6,
    "xtick.major.width": 0.6, "ytick.major.width": 0.6, "legend.frameon": False,
})



def main():
    with open(os.path.join(HERE, "retrospective_summary.json")) as f:
        s = json.load(f)

    with open(os.path.join(HERE, "context_summary.json")) as f:
        c = json.load(f)

    fig, (a, b) = plt.subplots(2, 1, figsize=(3.45, 4.3),
                               gridspec_kw={"height_ratios": [1.15, 1], "hspace": 0.55})

    # ---- (a) condition similarity vs pairwise identity -------------------
    labels = list(s["by_identity_bin"].keys())
    x = np.arange(len(labels))
    for off, series, color, name in [
        (-0.12, s["by_identity_bin"], BLUE, "All homolog pairs"),
        (0.0, c["same_composition_by_bin"], GREEN, "Same composition"),
        (0.12, s["uniprot"]["different_uniprot_by_bin"], ORANGE, "Different UniProt only"),
    ]:
        ys = [series[l]["jaccard"]["mean"] for l in labels]
        lo = [series[l]["jaccard"]["mean"] - series[l]["jaccard"]["ci95"][0] for l in labels]
        hi = [series[l]["jaccard"]["ci95"][1] - series[l]["jaccard"]["mean"] for l in labels]
        a.errorbar(x + off, ys, yerr=[lo, hi], color=color, lw=1.6, marker="o", ms=4,
                   mec="white", mew=0.8, capsize=0, elinewidth=1.0, label=name, zorder=3)
    base = s["baseline_unrelated"]["jaccard"]["mean"]
    a.axhline(base, color=INK2, lw=0.9, ls=(0, (3, 2)), zorder=1)
    a.text(len(labels) - 0.75, base - 0.012, "unrelated proteins", ha="right", va="top",
           color=INK2, fontsize=6.5)
    a.set_xticks(x, labels)
    a.set_xlabel("Pairwise sequence identity of the two entries")
    a.set_ylabel("Reagent-set Jaccard")
    a.set_ylim(0, 0.65)
    a.yaxis.grid(True, color=GRID, lw=0.5)
    a.set_axisbelow(True)
    a.legend(loc="upper left", fontsize=6.5, handlelength=1.6)
    a.set_title("(a)", loc="left", fontsize=7.5, color=INK)

    # ---- (b) same-UniProt pairs from different papers, by context ---------
    d = c["same_uniprot_different_paper"]
    rows = [("Same context", d["same_context"], GREEN),
            ("Apo vs. ligand-bound", d["same_context_apo_vs_ligand"], BLUE),
            ("Oligomeric state differs", d["only_olig_differs"], BLUE),
            ("Partner chains differ", d["only_partners_differ"], BLUE),
            ("Domain vs. longer construct", d["only_coverage_differs"], BLUE)]
    y = np.arange(len(rows))[::-1]
    for yy, (name, v, color) in zip(y, rows):
        m, (l, h) = v["jaccard"]["mean"], v["jaccard"]["ci95"]
        b.barh(yy, m, height=0.62, color=color, alpha=0.9 if color == GREEN else 0.55, zorder=2)
        b.errorbar(m, yy, xerr=[[m - l], [h - m]], color=INK2, lw=0.9, capsize=0, zorder=3)
        b.text(h + 0.012, yy, f"{m:.2f}  (n={v['jaccard']['n']:,})", va="center",
               fontsize=6, color=INK2)
    b.axvline(base, color=INK2, lw=0.9, ls=(0, (3, 2)), zorder=1)
    b.set_yticks(y, [r[0] for r in rows])
    b.set_xlim(0, 0.75)
    b.set_xlabel("Reagent-set Jaccard (same UniProt, different publications)")
    b.xaxis.grid(True, color=GRID, lw=0.5)
    b.set_axisbelow(True)
    b.spines["left"].set_visible(False)
    b.tick_params(axis="y", length=0)
    b.set_title("(b)", loc="left", fontsize=7.5, color=INK)

    fig.savefig(OUT + ".pdf", bbox_inches="tight")
    fig.savefig(OUT + ".png", dpi=300, bbox_inches="tight")
    print("wrote", OUT + ".pdf")


if __name__ == "__main__":
    main()
