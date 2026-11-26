"""Retrospective analysis on the held-out 10,000-entry panel: do homologs
crystallize in similar conditions, and would Crystal Explorer's homolog
conditions have recovered an entry's own (deposited) condition?

Inputs
  ../panel_10000_v2/panel_10000_raw.json   entries, target group, extractor output
  ../panel_10000_v2/targets.json           target sequences (one per group)
  entry_metadata.json                      from fetch_metadata.py (sequences,
                                           UniProt, pH, temperature, PubMed)
Outputs
  pairs.csv.gz                  every within-group pair with identity + similarities
  entries.csv                   per-entry derived fields
  retrospective_summary.json    all numbers quoted in the paper

Two analyses:

1. Pairwise (in the spirit of Abrahams & Newman 2019). Every pair of entries
   in the same target group (all hits of one query, i.e. homologs) gets a
   pairwise sequence identity (global alignment of the two matched chains,
   BLOSUM62, free end gaps, identities / length of the shorter chain), a
   same-UniProt flag and a same-publication flag, and three condition
   similarities: Jaccard of reagent sets, shared precipitant, |delta pH|.
   Random pairs of entries from different target groups (different
   30%-identity clusters, i.e. unrelated proteins) are the baseline.

2. Leave-one-out retrieval. Each entry is treated as a "new" protein: its
   own record is hidden and the k most similar remaining entries of its
   group (ranked by identity to it) play the role of Crystal Explorer's
   suggestions. Three scenarios remove progressively more information:
   all homologs; homologs from other publications; homologs with a
   different UniProt accession (a genuinely different protein). Baselines
   use the same k: k random entries from unrelated groups, and a
   "generic screen" made of the globally most frequent reagents.

Uncertainty: 95% CIs by cluster bootstrap over target groups (pairs and
queries within one group are not independent), 1,000 resamples.
"""
import gzip
import json
import os
import re
from collections import Counter, defaultdict
from itertools import combinations

import numpy as np
import pandas as pd
from Bio import Align
from Bio.Align import substitution_matrices

HERE = os.path.dirname(os.path.abspath(__file__))
PANEL_DIR = os.path.join(HERE, "..", "panel_10000_v2")
SEED = 2030
N_BOOT = 1000
K = 5                 # number of homolog conditions shown to the user (top-k)
N_RANDOM_PAIRS = 200_000
SALT_PRECIP_M = 0.8   # a salt counts as precipitant at >= this molarity

rng = np.random.default_rng(SEED)

# --------------------------------------------------------------------------
# sequence identity
# --------------------------------------------------------------------------
ALIGNER = Align.PairwiseAligner()
ALIGNER.substitution_matrix = substitution_matrices.load("BLOSUM62")
ALIGNER.open_gap_score = -10
ALIGNER.extend_gap_score = -0.5
ALIGNER.end_gap_score = 0
_ALPHA = set(ALIGNER.substitution_matrix.alphabet)


def clean(seq):
    return "".join(c if c in _ALPHA else "X" for c in (seq or "").upper())


def identity(a, b):
    """Identities / length of the shorter chain, in percent."""
    if not a or not b:
        return np.nan
    aln = ALIGNER.align(a, b)[0]
    return 100.0 * aln.counts().identities / min(len(a), len(b))


# --------------------------------------------------------------------------
# reagents
# --------------------------------------------------------------------------
ITEM_RE = re.compile(r"'([^']+)'\s*\(([^)]*)\)")
MERGE = {  # spelling / counter-ion variants of the same reagent
    "tris-hcl": "tris", "tris hcl": "tris", "tris base": "tris",
    "sodium hepes": "hepes", "hepes-naoh": "hepes",
    "trisodium citrate": "sodium citrate", "citrate buffer": "sodium citrate",
    "sodium/potassium phosphate": "sodium/potassium phosphate",
    "2-propanol": "isopropanol",
}
SALT_WORDS = ("sulfate", "sulphate", "chloride", "phosphate", "citrate", "malonate",
              "formate", "acetate", "tartrate", "nitrate", "thiocyanate", "bromide",
              "iodide", "fluoride", "succinate", "hpo4", "po4", "so4")
ORGANIC_PRECIP = ("mpd", "isopropanol", "ethanol", "hexanediol", "jeffamine",
                  "pentaerythritol", "tert-butanol", "methanol", "butanediol")


def norm_name(n):
    n = n.strip().lower()
    return MERGE.get(n, n)


def peg_class(n):
    if "peg" not in n and "polyethylene glycol" not in n:
        return None
    nums = [int(x) for x in re.findall(r"\d+", n.replace(",", ""))]
    nums = [x for x in nums if x >= 100]
    if not nums:
        return "PEG (unspecified MW)"
    mw = nums[0]
    if mw < 1000:
        return "PEG <1k"
    if mw < 6000:
        return "PEG 1-5k"
    return "PEG >=6k"


def molarity(conc):
    m = re.match(r"\s*([\d.]+)(?:\s*-\s*([\d.]+))?\s*(mM|M)\b", conc or "")
    if not m:
        return None
    lo = float(m.group(1))
    hi = float(m.group(2)) if m.group(2) else lo
    v = (lo + hi) / 2
    return v / 1000 if m.group(3) == "mM" else v


def parse_reagents(s):
    reagents, precip = set(), set()
    for name, conc in ITEM_RE.findall(s or ""):
        n = norm_name(name)
        reagents.add(n)
        pc = peg_class(n)
        if pc:
            precip.add(pc)
        elif any(w in n for w in ORGANIC_PRECIP) or n == "tacsimate":
            precip.add(n)
        elif any(w in n for w in SALT_WORDS):
            m = molarity(conc)
            if m is not None and m >= SALT_PRECIP_M:
                precip.add(n)
    return frozenset(reagents), frozenset(precip)


def jaccard(a, b):
    if not a and not b:
        return np.nan
    return len(a & b) / len(a | b)


# --------------------------------------------------------------------------
# cluster bootstrap helpers
# --------------------------------------------------------------------------
def boot_groups(groups):
    """List of index arrays, one per bootstrap resample of target groups."""
    uniq = np.unique(groups)
    idx_by_g = {g: np.flatnonzero(groups == g) for g in uniq}
    out = []
    for _ in range(N_BOOT):
        pick = rng.choice(uniq, size=len(uniq), replace=True)
        out.append(np.concatenate([idx_by_g[g] for g in pick]))
    return out


def mean_ci(x, groups, boots=None):
    x = np.asarray(x, float)
    ok = ~np.isnan(x)
    if ok.sum() == 0:
        return {"n": 0}
    x, groups = x[ok], np.asarray(groups)[ok]
    boots = boots if boots is not None else boot_groups(groups)
    bs = [np.nanmean(x[b]) for b in boots]
    return {"n": int(ok.sum()), "mean": round(float(x.mean()), 4),
            "ci95": [round(float(np.percentile(bs, 2.5)), 4), round(float(np.percentile(bs, 97.5)), 4)]}


def median_ci(x, groups):
    x = np.asarray(x, float)
    ok = ~np.isnan(x)
    x, groups = x[ok], np.asarray(groups)[ok]
    bs = [np.median(x[b]) for b in boot_groups(groups)]
    return {"n": int(ok.sum()), "median": round(float(np.median(x)), 4),
            "ci95": [round(float(np.percentile(bs, 2.5)), 4), round(float(np.percentile(bs, 97.5)), 4)]}


def rankdata(x):
    """Average ranks (ties share the mean rank)."""
    x = np.asarray(x)
    order = np.argsort(x, kind="mergesort")
    ranks = np.empty(len(x), float)
    sx = x[order]
    bounds = np.flatnonzero(np.diff(sx)) + 1
    starts = np.concatenate([[0], bounds])
    ends = np.concatenate([bounds, [len(x)]])
    for s, e in zip(starts, ends):
        ranks[order[s:e]] = (s + e - 1) / 2.0 + 1
    return ranks


def spearman(x, y):
    rx, ry = rankdata(x), rankdata(y)
    return float(np.corrcoef(rx, ry)[0, 1])


def spearman_ci(x, y, groups):
    x, y, groups = map(np.asarray, (x, y, groups))
    ok = ~(np.isnan(x) | np.isnan(y))
    x, y, groups = x[ok], y[ok], groups[ok]
    rho = spearman(x, y)
    bs = [spearman(x[b], y[b]) for b in boot_groups(groups)]
    # one-sided permutation p-value: shuffle y within target groups
    null = []
    gidx = [np.flatnonzero(groups == g) for g in np.unique(groups)]
    for _ in range(200):
        yp = y.copy()
        for ix in gidx:
            yp[ix] = y[rng.permutation(ix)]
        null.append(spearman(x, yp))
    p = (1 + sum(abs(n) >= abs(rho) for n in null)) / (1 + len(null))
    return {"n_pairs": int(ok.sum()), "rho": round(rho, 4),
            "ci95": [round(float(np.percentile(bs, 2.5)), 4), round(float(np.percentile(bs, 97.5)), 4)],
            "perm_p_within_group": round(p, 4)}


# --------------------------------------------------------------------------
def build_entries():
    """One row per panel entry: matched chain, UniProt, PubMed, pH, deposit
    date, reagent and precipitant sets."""
    with open(os.path.join(PANEL_DIR, "panel_10000_raw.json")) as f:
        panel = json.load(f)
    with open(os.path.join(PANEL_DIR, "targets.json")) as f:
        targets = json.load(f)
    with open(os.path.join(HERE, "entry_metadata.json")) as f:
        meta = json.load(f)
    rows = []
    for e in panel:
        pid = e["pdb_id"].upper()
        m = meta[pid]
        tseq = clean(targets[e["protein"]]["sequence"])
        best = None
        for ent in m["entities"]:
            if not (ent["type"] or "").startswith("polypeptide"):
                continue
            s = clean(ent["seq"])
            idt = identity(tseq, s)
            if best is None or idt > best[0]:
                best = (idt, s, ent)
        reag, prec = parse_reagents(e["extracted_compound_string"])
        rows.append({
            "pdb_id": pid, "group": e["protein"],
            "seq": best[1] if best else "", "id_to_target": best[0] if best else np.nan,
            "uniprot": frozenset(best[2]["uniprot"]) if best else frozenset(),
            "pubmed": m["pubmed"], "pH": m["pH"], "temp": m["temp"],
            "deposit_date": (m.get("deposit_date") or "")[:10],
            "reagents": reag, "precip": prec,
        })
    return pd.DataFrame(rows)


def main():
    E = build_entries()
    E.to_csv(os.path.join(HERE, "entries.csv"), index=False,
             columns=["pdb_id", "group", "id_to_target", "pubmed", "pH", "temp"])
    print(f"{len(E)} entries, {E.group.nunique()} groups, "
          f"{(E.reagents.map(len) > 0).sum()} with reagents, "
          f"{(E.precip.map(len) > 0).sum()} with an identified precipitant, "
          f"{E.uniprot.map(len).gt(0).sum()} with UniProt")

    # ---- within-group pairs ----------------------------------------------
    pairs = []
    gid = {g: k for k, g in enumerate(sorted(E.group.unique()))}
    for g, sub in E.groupby("group"):
        idx = sub.index.tolist()
        for i, j in combinations(idx, 2):
            a, b = E.loc[i], E.loc[j]
            pairs.append((i, j, gid[g], identity(a.seq, b.seq)))
    P = pd.DataFrame(pairs, columns=["i", "j", "g", "identity"])

    def pair_feats(P):
        A, B = E.loc[P.i].reset_index(drop=True), E.loc[P.j].reset_index(drop=True)
        P = P.reset_index(drop=True)
        P["jaccard"] = [jaccard(a, b) if (a and b) else np.nan for a, b in zip(A.reagents, B.reagents)]
        P["shared_precip"] = [float(bool(a & b)) if (a and b) else np.nan for a, b in zip(A.precip, B.precip)]
        P["dpH"] = (A.pH.astype(float) - B.pH.astype(float)).abs()
        P["same_uniprot"] = [bool(a & b) for a, b in zip(A.uniprot, B.uniprot)]
        P["uniprot_known"] = [bool(a) and bool(b) for a, b in zip(A.uniprot, B.uniprot)]
        P["same_paper"] = [(a is not None and a == b) for a, b in zip(A.pubmed, B.pubmed)]
        return P

    P = pair_feats(P)
    with gzip.open(os.path.join(HERE, "pairs.csv.gz"), "wt") as f:
        P.assign(pdb_i=E.pdb_id.values[P.i], pdb_j=E.pdb_id.values[P.j]).to_csv(f, index=False)
    print(f"{len(P)} within-group pairs")

    # random cross-group pairs = unrelated proteins
    gi = E.group.map(gid).values
    ri = rng.integers(0, len(E), N_RANDOM_PAIRS)
    rj = rng.integers(0, len(E), N_RANDOM_PAIRS)
    keep = gi[ri] != gi[rj]
    R = pair_feats(pd.DataFrame({"i": ri[keep], "j": rj[keep], "g": gi[ri[keep]], "identity": np.nan}))
    R["g"] = np.arange(len(R))  # independent pairs -> ordinary bootstrap

    summary = {"params": {"seed": SEED, "n_boot": N_BOOT, "k": K, "salt_precip_M": SALT_PRECIP_M},
               "entries": int(len(E)), "groups": int(E.group.nunique()),
               "within_group_pairs": int(len(P)), "random_pairs": int(len(R))}

    metrics = ["jaccard", "shared_precip", "dpH"]
    summary["baseline_unrelated"] = {m: mean_ci(R[m], R.g) for m in metrics}

    bins = [(0, 50), (50, 70), (70, 90), (90, 98), (98, 100.01)]
    labels = ["<50%", "50-70%", "70-90%", "90-98%", ">=98%"]
    P["bin"] = pd.cut(P.identity, [b[0] for b in bins] + [bins[-1][1]], right=False, labels=labels)

    def strata(df):
        out = {}
        for m in metrics:
            out[m] = mean_ci(df[m], df.g)
        return out

    summary["all_homolog_pairs"] = strata(P)
    summary["by_identity_bin"] = {l: strata(P[P.bin == l]) for l in labels}
    U = P[P.uniprot_known]
    summary["uniprot"] = {
        "same_uniprot": strata(U[U.same_uniprot]),
        "different_uniprot": strata(U[~U.same_uniprot]),
        "same_uniprot_different_paper": strata(U[U.same_uniprot & ~U.same_paper]),
        "same_paper": strata(P[P.same_paper]),
        # does UniProt identity matter beyond sequence identity? compare at >=90%
        "ge90_same_uniprot": strata(U[(U.identity >= 90) & U.same_uniprot]),
        "ge90_different_uniprot": strata(U[(U.identity >= 90) & ~U.same_uniprot]),
        "different_uniprot_by_bin": {l: strata(U[(~U.same_uniprot) & (U.bin == l)]) for l in labels},
        "same_uniprot_fraction_by_bin": {l: round(float(U[U.bin == l].same_uniprot.mean()), 4)
                                         for l in labels if (U.bin == l).any()},
    }
    summary["bin_counts"] = {l: int((P.bin == l).sum()) for l in labels}

    summary["spearman_identity_vs"] = {
        "all_pairs": {m: spearman_ci(P.identity, P[m] if m != "dpH" else -P[m], P.g) for m in metrics},
        "different_uniprot": {m: spearman_ci(U[~U.same_uniprot].identity,
                                             U[~U.same_uniprot][m] if m != "dpH" else -U[~U.same_uniprot][m],
                                             U[~U.same_uniprot].g) for m in metrics},
        "different_paper": {m: spearman_ci(P[~P.same_paper].identity,
                                           P[~P.same_paper][m] if m != "dpH" else -P[~P.same_paper][m],
                                           P[~P.same_paper].g) for m in metrics},
        "note": "dpH is negated so a positive rho always means 'more identical -> more similar conditions'",
    }

    # ---- leave-one-out retrieval -----------------------------------------
    pid_idx = {}
    for r in P.itertuples():
        pid_idx[(r.i, r.j)] = r.identity
        pid_idx[(r.j, r.i)] = r.identity
    members = defaultdict(list)
    for i, g in enumerate(E.group):
        members[g].append(i)

    freq = Counter(r for s in E.reagents for r in s)
    freq_p = Counter(r for s in E.precip for r in s)
    ranked_reag = [r for r, _ in freq.most_common()]
    ranked_prec = [r for r, _ in freq_p.most_common()]
    global_pH = float(np.nanmedian(E.pH.astype(float)))
    has_reag = E.reagents.map(len).values > 0
    all_with_reag = np.flatnonzero(has_reag)
    unrelated_pool = {g: all_with_reag[E.group.values[all_with_reag] != g] for g in members}

    def suggest_metrics(q, cand):
        """Metrics for query q given a list of suggested entry indices."""
        Rq, Pq = E.reagents[q], E.precip[q]
        U_r = frozenset().union(*[E.reagents[c] for c in cand])
        U_p = frozenset().union(*[E.precip[c] for c in cand])
        top1 = cand[0]
        pHs = np.array([E.pH[c] for c in cand if pd.notna(E.pH[c])], float)
        pq = E.pH[q] if pd.notna(E.pH[q]) else None
        return {
            "top1_jaccard": jaccard(Rq, E.reagents[top1]),
            "topk_best_jaccard": max(jaccard(Rq, E.reagents[c]) for c in cand),
            "topk_precip_recovered": float(bool(Pq & U_p)) if Pq else np.nan,
            "topk_reagent_recall": len(Rq & U_r) / len(Rq),
            "suggested_set_size": len(U_r),
            "n_precip_suggested": len(U_p),
            "abs_dpH_median": abs(pq - np.median(pHs)) if (pq is not None and len(pHs)) else np.nan,
            "pH_in_range": float(pHs.min() - 0.25 <= pq <= pHs.max() + 0.25) if (pq is not None and len(pHs)) else np.nan,
        }

    def generic_metrics(q, n_reag, n_prec):
        Rq, Pq = E.reagents[q], E.precip[q]
        U_r = set(ranked_reag[:n_reag])
        U_p = set(ranked_prec[:max(n_prec, 1)])
        pq = E.pH[q] if pd.notna(E.pH[q]) else None
        return {
            "topk_reagent_recall": len(Rq & U_r) / len(Rq),
            "topk_precip_recovered": float(bool(Pq & U_p)) if Pq else np.nan,
            "abs_dpH_median": abs(pq - global_pH) if pq is not None else np.nan,
        }

    scenarios = {
        "all_homologs": lambda q, c: True,
        "other_papers": lambda q, c: not (E.pubmed[q] is not None and E.pubmed[q] == E.pubmed[c]),
        "different_uniprot": lambda q, c: bool(E.uniprot[q]) and bool(E.uniprot[c])
                                          and not (E.uniprot[q] & E.uniprot[c]),
    }
    loo = {}
    for sname, allowed in scenarios.items():
        recs = []
        for q in range(len(E)):
            if not has_reag[q]:
                continue
            g = E.group[q]
            cand = [c for c in members[g] if c != q and has_reag[c] and allowed(q, c)]
            if not cand:
                continue
            cand.sort(key=lambda c: -pid_idx[(q, c)])
            cand = cand[:K]
            hm = suggest_metrics(q, cand)
            # random baseline: same number of entries from unrelated groups, 5 draws averaged
            rnd = []
            pool = unrelated_pool[g]
            for _ in range(5):
                rc = list(rng.choice(pool, size=len(cand), replace=False))
                rnd.append(suggest_metrics(q, rc))
            rm = {k: np.nanmean([d[k] for d in rnd]) if not all(np.isnan(d[k]) for d in rnd) else np.nan
                  for k in rnd[0]}
            gm = generic_metrics(q, hm["suggested_set_size"], hm["n_precip_suggested"])
            recs.append({"q": q, "g": g, "k": len(cand),
                         **{f"hom_{k}": v for k, v in hm.items()},
                         **{f"rnd_{k}": v for k, v in rm.items()},
                         **{f"gen_{k}": v for k, v in gm.items()}})
        D = pd.DataFrame(recs)
        D.to_csv(os.path.join(HERE, f"loo_{sname}.csv"), index=False)
        boots = boot_groups(D.g.values)
        res = {"queries": int(len(D)), "groups": int(D.g.nunique()),
               "mean_k": round(float(D.k.mean()), 2)}
        for metric in ["top1_jaccard", "topk_best_jaccard", "topk_precip_recovered",
                       "topk_reagent_recall", "abs_dpH_median", "pH_in_range"]:
            res[metric] = {"homologs": mean_ci(D[f"hom_{metric}"], D.g, None),
                           "random_unrelated": mean_ci(D[f"rnd_{metric}"], D.g, None)}
            if f"gen_{metric}" in D:
                res[metric]["generic_frequent"] = mean_ci(D[f"gen_{metric}"], D.g, None)
        # median |dpH| is more robust than the mean
        res["abs_dpH_median"]["homologs_median"] = median_ci(D["hom_abs_dpH_median"], D.g)
        res["abs_dpH_median"]["random_median"] = median_ci(D["rnd_abs_dpH_median"], D.g)
        res["abs_dpH_median"]["generic_median"] = median_ci(D["gen_abs_dpH_median"], D.g)
        # fraction of queries where homolog top-1 equals the query's reagent set exactly
        res["top1_exact_reagent_match"] = mean_ci((D.hom_top1_jaccard == 1).astype(float), D.g)
        loo[sname] = res
        print(sname, res["queries"], "queries;",
              "precip recovered", res["topk_precip_recovered"]["homologs"]["mean"],
              "vs random", res["topk_precip_recovered"]["random_unrelated"]["mean"],
              "vs generic", res["topk_precip_recovered"]["generic_frequent"]["mean"])
    summary["leave_one_out"] = loo
    summary["global_median_pH"] = global_pH

    with open(os.path.join(HERE, "retrospective_summary.json"), "w") as f:
        json.dump(summary, f, indent=2, default=lambda o: None)
    print(json.dumps({k: summary[k] for k in ["baseline_unrelated", "by_identity_bin", "spearman_identity_vs"]},
                     indent=1))


if __name__ == "__main__":
    main()
