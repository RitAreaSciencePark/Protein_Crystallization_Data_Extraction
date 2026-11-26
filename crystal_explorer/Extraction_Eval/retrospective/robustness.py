"""Two robustness checks on the retrospective analysis (paper Section VI-D).

1. Per-target weighting. In the main analysis every pair / query counts
   once, so targets with long series of near-identical structures weigh
   more. Here every target counts once: metrics are averaged within a
   target first, then across targets (cluster bootstrap CIs over targets).

2. Time split. Entries deposited on or after CUTOFF are queries; only
   entries deposited before CUTOFF may be suggested, and the generic
   frequency list is built from those earlier entries only. This mimics
   real use (nobody sees future structures) and removes the small leakage
   of the main analysis, where the generic list included the query.

Needs analyze_retrospective.py to have been run (pairs.csv.gz, loo_*.csv).

    python robustness.py   # -> robustness_summary.json
"""
import gzip
import json
import os
from collections import Counter, defaultdict

import numpy as np
import pandas as pd

from analyze_retrospective import K, N_BOOT, build_entries, jaccard

HERE = os.path.dirname(os.path.abspath(__file__))
CUTOFF = "2020-01-01"
SEED = 2031
rng = np.random.default_rng(SEED)
SCENARIOS = ["all_homologs", "other_papers", "different_uniprot"]


def macro(df, col, group="g"):
    """Mean of per-target means, with a bootstrap CI over targets."""
    d = df[[group, col]].dropna()
    if d.empty:
        return {"targets": 0}
    per = d.groupby(group)[col].mean().values
    bs = [rng.choice(per, len(per), replace=True).mean() for _ in range(N_BOOT)]
    return {"targets": int(len(per)), "mean": round(float(per.mean()), 4),
            "ci95": [round(float(np.percentile(bs, 2.5)), 4), round(float(np.percentile(bs, 97.5)), 4)]}


def micro(df, col, group="g"):
    d = df[[group, col]].dropna()
    groups = d[group].values
    uniq = np.unique(groups)
    idx = {g: np.flatnonzero(groups == g) for g in uniq}
    x = d[col].values.astype(float)
    bs = []
    for _ in range(N_BOOT):
        pick = np.concatenate([idx[g] for g in rng.choice(uniq, len(uniq), replace=True)])
        bs.append(x[pick].mean())
    return {"n": int(len(x)), "targets": int(len(uniq)), "mean": round(float(x.mean()), 4),
            "ci95": [round(float(np.percentile(bs, 2.5)), 4), round(float(np.percentile(bs, 97.5)), 4)]}


def median_of(df, col):
    d = df[col].dropna()
    return round(float(d.median()), 3) if len(d) else None


def main():
    summary = {"cutoff": CUTOFF, "seed": SEED, "k": K}

    # ------------------------------------------------------------------
    # 1. per-target weighting
    # ------------------------------------------------------------------
    with gzip.open(os.path.join(HERE, "pairs.csv.gz"), "rt") as f:
        P = pd.read_csv(f)
    bins = [0, 50, 70, 90, 98, 100.01]
    labels = ["<50%", "50-70%", "70-90%", "90-98%", ">=98%"]
    P["bin"] = pd.cut(P.identity, bins, right=False, labels=labels)
    U = P[P.uniprot_known]
    pw = {"by_identity_bin": {l: macro(P[P.bin == l], "jaccard") for l in labels},
          "same_uniprot": macro(U[U.same_uniprot], "jaccard"),
          "different_uniprot": macro(U[~U.same_uniprot], "jaccard"),
          "same_uniprot_different_paper": macro(U[U.same_uniprot & ~U.same_paper], "jaccard")}
    # how much of the pair mass sits in a few targets
    per_target = P.groupby("g").size().sort_values(ascending=False)
    pw["share_of_pairs_in_top10pct_targets"] = round(
        float(per_target.iloc[:max(1, len(per_target) // 10)].sum() / per_target.sum()), 3)
    summary["per_target_weighting"] = {"pairs_jaccard": pw}

    loo = {}
    for s in SCENARIOS:
        D = pd.read_csv(os.path.join(HERE, f"loo_{s}.csv"))
        loo[s] = {m: {src: macro(D, f"{src}_{m}") for src in ("hom", "gen", "rnd")}
                  for m in ("topk_precip_recovered", "topk_reagent_recall")}
    summary["per_target_weighting"]["leave_one_out"] = loo

    # ------------------------------------------------------------------
    # 2. time split
    # ------------------------------------------------------------------
    E = build_entries()
    date = E.deposit_date.fillna("")
    is_train = (date < CUTOFF) & (date != "")
    is_test = date >= CUTOFF
    has_reag = E.reagents.map(len).values > 0
    summary["time_split_sizes"] = {"train_entries": int(is_train.sum()), "test_entries": int(is_test.sum()),
                                   "test_with_reagents": int((is_test & has_reag).sum())}

    ident = {}
    for r in P.itertuples():
        ident[(r.i, r.j)] = r.identity
        ident[(r.j, r.i)] = r.identity
    train_by_group = defaultdict(list)
    for i in np.flatnonzero(is_train.values & has_reag):
        train_by_group[E.group[i]].append(i)
    train_all = np.flatnonzero(is_train.values & has_reag)
    train_groups = E.group.values[train_all]

    freq = Counter(r for i in train_all for r in E.reagents[i])
    freq_p = Counter(p for i in train_all for p in E.precip[i])
    ranked_r = [r for r, _ in freq.most_common()]
    ranked_p = [p for p, _ in freq_p.most_common()]
    train_pH = float(np.nanmedian(E.pH.values[train_all].astype(float)))

    def ph(i):
        v = E.pH[i]
        return None if pd.isna(v) else float(v)

    def metrics(q, cand):
        Rq, Pq = E.reagents[q], E.precip[q]
        Ur = frozenset().union(*[E.reagents[c] for c in cand])
        Up = frozenset().union(*[E.precip[c] for c in cand])
        pHs = [ph(c) for c in cand if ph(c) is not None]
        pq = ph(q)
        return {"precip": float(bool(Pq & Up)) if Pq else np.nan,
                "recall": len(Rq & Ur) / len(Rq),
                "dpH": abs(pq - float(np.median(pHs))) if (pq is not None and pHs) else np.nan,
                "top1_jaccard": jaccard(Rq, E.reagents[cand[0]]),
                "n_r": len(Ur), "n_p": len(Up)}

    allowed = {
        "all_homologs": lambda q, c: True,
        "other_papers": lambda q, c: not (E.pubmed[q] is not None and not pd.isna(E.pubmed[q])
                                          and E.pubmed[q] == E.pubmed[c]),
        "different_uniprot": lambda q, c: bool(E.uniprot[q]) and bool(E.uniprot[c])
                                          and not (E.uniprot[q] & E.uniprot[c]),
    }
    ts = {}
    for s in SCENARIOS:
        recs = []
        for q in np.flatnonzero(is_test.values & has_reag):
            g = E.group[q]
            cand = [c for c in train_by_group[g] if allowed[s](q, c)]
            if not cand:
                continue
            cand.sort(key=lambda c: -ident.get((q, c), -1))
            cand = cand[:K]
            h = metrics(q, cand)
            pool = train_all[train_groups != g]
            rnd = [metrics(q, list(rng.choice(pool, len(cand), replace=False))) for _ in range(5)]
            Pq, Rq = E.precip[q], E.reagents[q]
            gen_p = float(bool(Pq & set(ranked_p[:max(h["n_p"], 1)]))) if Pq else np.nan
            gen_r = len(Rq & set(ranked_r[:h["n_r"]])) / len(Rq)
            pq = ph(q)
            recs.append({"g": g, "k": len(cand),
                         "top1_same_uniprot": bool(E.uniprot[q] & E.uniprot[cand[0]]),
                         "hom_precip": h["precip"], "hom_recall": h["recall"], "hom_dpH": h["dpH"],
                         "hom_top1_jaccard": h["top1_jaccard"],
                         "rnd_precip": np.nanmean([r["precip"] for r in rnd]) if Pq else np.nan,
                         "rnd_recall": np.mean([r["recall"] for r in rnd]),
                         "rnd_dpH": np.nanmean([r["dpH"] for r in rnd]) if pq is not None else np.nan,
                         "gen_precip": gen_p, "gen_recall": gen_r,
                         "gen_dpH": abs(pq - train_pH) if pq is not None else np.nan})
        D = pd.DataFrame(recs)
        D.to_csv(os.path.join(HERE, f"timesplit_{s}.csv"), index=False)
        res = {"queries": int(len(D)), "targets": int(D.g.nunique()), "mean_k": round(float(D.k.mean()), 2)}
        for m in ("precip", "recall"):
            res[m] = {src: micro(D, f"{src}_{m}") for src in ("hom", "gen", "rnd")}
            res[m + "_per_target"] = {src: macro(D, f"{src}_{m}") for src in ("hom", "gen", "rnd")}
        res["median_abs_dpH"] = {src: median_of(D, f"{src}_dpH") for src in ("hom", "gen", "rnd")}
        res["top1_exact_reagent_match"] = round(float((D.hom_top1_jaccard == 1).mean()), 4)
        res["share_top1_same_uniprot"] = round(float(D.top1_same_uniprot.mean()), 4)
        res["by_top1_same_uniprot"] = {
            str(flag): {"queries": int(len(sub)),
                        **{f"{src}_{m}": micro(sub, f"{src}_{m}") for src in ("hom", "gen") for m in ("precip", "recall")}}
            for flag, sub in D.groupby("top1_same_uniprot")}
        for a, b, m in (("hom", "gen", "precip"), ("hom", "gen", "recall"), ("gen", "hom", "dpH")):
            d = D[["g", f"{a}_{m}", f"{b}_{m}"]].dropna()
            d = d.assign(x=d[f"{a}_{m}"] - d[f"{b}_{m}"])
            res[f"paired_diff_{a}_minus_{b}_{m}"] = micro(d, "x")
        ts[s] = res
        print(f"time split {s}: {len(D)} queries / {D.g.nunique()} targets; precip "
              f"hom {res['precip']['hom']['mean']} gen {res['precip']['gen']['mean']} "
              f"rnd {res['precip']['rnd']['mean']}; recall hom {res['recall']['hom']['mean']} "
              f"gen {res['recall']['gen']['mean']}")
    summary["time_split"] = ts

    with open(os.path.join(HERE, "robustness_summary.json"), "w") as f:
        json.dump(summary, f, indent=2)
    print(json.dumps(summary["per_target_weighting"], indent=1))


if __name__ == "__main__":
    main()
