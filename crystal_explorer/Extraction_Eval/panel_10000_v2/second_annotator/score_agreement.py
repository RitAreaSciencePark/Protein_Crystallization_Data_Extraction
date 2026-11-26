"""Inter-annotator agreement and extractor scores against a blinded second
annotator (paper Section VI-B).

Inputs (all in this directory)
  annotator_A_structured.csv              A's ground truth, reviewed by A (no "[FN?]" left)
  for_annotator_B/annotation_sheet.csv    filled in by annotator B (returned file)
  ../audit_1000_raw.json                  extractor output at audit time

Two steps:

  python score_agreement.py prepare
      Matches names automatically (same chemical after normalization) and
      writes adjudication.csv with every item that could not be matched
      automatically, for A vs B, extractor vs B, and each LLM run vs B
      (every ../llm_comparison/llm_runs/<model>/predictions.jsonl that covers
      all 200 entries). Re-running keeps any decisions already typed in, so
      a new LLM run can be added later and only its new rows need filling.

  -> Fill in adjudication.csv, ideally A and B together:
       same_as    if the item is the same reagent as one listed in
                  `other_side_unmatched` under the paper's scoring rules
                  (including a right-family name, e.g. generic "PEG" for
                  "PEG 8000"), copy that name here exactly; else leave blank.
       consensus  for A_vs_B rows left unmatched only: "reagent" if the item
                  should be in the ground truth, "not_reagent" if not.

  python score_agreement.py score
      -> agreement_summary.json, disagreements.csv

Metrics (95% CIs by bootstrap over entries, 2,000 resamples)
  A vs B   positive specific agreement (PSA = 2*matches / (|A|+|B|), the
           mention-level F1 between annotators; kappa is undefined here
           because there are no true negatives) and entry-level exact
           agreement.
  extractor vs A           the original audit counts on these 200 entries
  extractor vs B           precision / recall / F1 / entry exact match
  extractor vs consensus   same, against the adjudicated ground truth
  llm:<model> vs B / vs consensus   the same, for each finished LLM run
"""
import csv
import glob
import json
import os
import re
import sys

import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))
PANEL = os.path.dirname(HERE)
ROOT = os.path.abspath(os.path.join(HERE, "..", "..", "..", ".."))
DICT_PATH = os.path.join(ROOT, "Input", "compound_dictionary.json")
ADJ = os.path.join(HERE, "adjudication.csv")
LLM_RUNS = os.path.join(PANEL, "llm_comparison", "llm_runs")
BASELINES = os.path.join(PANEL, "baselines")
ITEM_RE = re.compile(r"'([^']+)'\s*\(")
N_BOOT = 2000
rng = np.random.default_rng(2033)

with open(DICT_PATH) as f:
    ALIAS = json.load(f)

BUFFERS = {"tris", "hepes", "mes", "mops", "pipes", "caps", "ches", "bicine", "tricine",
           "bistris", "bistrispropane", "imidazole", "cacodylate", "citrate", "acetate",
           "glycine", "taps", "epps", "heppso", "mopso", "ada", "aces", "bes", "tes"}
EXTRA = {  # same chemical, different spelling the dictionary does not unify
    "dithiothreitol": "dtt", "betamercaptoethanol": "2mercaptoethanol", "bme": "2mercaptoethanol",
    "mercaptoethanol": "2mercaptoethanol", "zinciiacetate": "zincacetate",
    "hepesna": "hepes", "sodiumhepes": "hepes", "nahepes": "hepes",
    "trishydrochloride": "tris", "trisbase": "tris",
}


def squash(s):
    s = s.lower().replace("sulphate", "sulfate").replace("polyethylene glycol", "peg")
    s = s.replace("monomethyl ether", "mme").replace("monomethylether", "mme")
    return re.sub(r"[^a-z0-9]", "", s)


def key(name):
    """Normalized identity of a reagent name, for automatic matching only."""
    k = squash(re.sub(r"\[fn\?\]", "", name, flags=re.I))
    if k in ALIAS:
        k = squash(ALIAS[k])
    k = EXTRA.get(k, k)
    for suf in ("hcl", "naoh", "koh", "na"):
        if k.endswith(suf) and k[: -len(suf)] in BUFFERS:
            k = k[: -len(suf)]
    return EXTRA.get(k, k)


def split_list(s):
    s = (s or "").strip()
    if not s or s.upper() == "NONE":
        return []
    out = []
    for x in s.split(";"):
        x = x.strip()
        if x and x not in out:
            out.append(x)
    return out


def load():
    with open(os.path.join(HERE, "annotator_A_structured.csv"), encoding="utf-8-sig") as f:
        A = {r["pdb_id"]: r for r in csv.DictReader(f)}
    with open(os.path.join(HERE, "for_annotator_B", "annotation_sheet.csv"), encoding="utf-8-sig") as f:
        B = {r["pdb_id"]: r for r in csv.DictReader(f)}
    with open(os.path.join(PANEL, "audit_1000_raw.json")) as f:
        X = {e["pdb_id"]: e for e in json.load(f)}
    left = [p for p, r in A.items() if "[FN?]" in r["A_reagents"]]
    if left:
        sys.exit(f"annotator_A_structured.csv still has unreviewed [FN?] items: {left}")
    empty = [p for p, r in B.items() if not r["reagents"].strip()]
    if empty:
        sys.exit(f"annotator B left {len(empty)} entries blank (write NONE if no reagent): {empty[:10]}")
    order = sorted(A, key=lambda p: int(A[p]["item"]))
    data = {}
    for p in order:
        ext = list(dict.fromkeys(ITEM_RE.findall(X[p]["extracted_compound_string"] or "")))
        data[p] = {"item": int(A[p]["item"]), "details": X[p]["pdbx_details"],
                   "A": split_list(A[p]["A_reagents"]), "B": split_list(B[p]["reagents"]),
                   "systems": {"extractor": ext},
                   "orig": (int(A[p]["tp"]), int(A[p]["fp"]), int(A[p]["fn"]))}
    for name, preds in load_llm_runs(set(data)).items():
        for p in data:
            data[p]["systems"][name] = preds[p]
    return data


def load_llm_runs(ids):
    """Every other system's predictions that cover all sampled entries without
    error: LLM runs (../llm_comparison/llm_runs/<model>/, from llm_extract.py)
    as "llm:<model>", and baselines (../baselines/<name>/, e.g. the Lynch et
    al. 2020 parser from run_lynch2020.py) as "baseline:<name>".
    Returns {system: {pdb_id: [reagents]}}."""
    runs = {}
    paths = [("llm:", p) for p in sorted(glob.glob(os.path.join(LLM_RUNS, "*", "predictions.jsonl")))] + \
            [("baseline:", p) for p in sorted(glob.glob(os.path.join(BASELINES, "*", "predictions.jsonl")))]
    for prefix, path in paths:
        preds = {}
        with open(path) as f:
            for line in f:
                r = json.loads(line)
                if "error" not in r:
                    preds[r["pdb_id"]] = r["reagents"]
        name = prefix + os.path.basename(os.path.dirname(path))
        missing = ids - set(preds)
        if missing:
            print(f"skipping {name}: {len(missing)} of {len(ids)} entries not (successfully) extracted")
            continue
        runs[name] = preds
    return runs


def auto_match(left, right):
    """Greedy one-to-one match on normalized key. Returns (pairs, unmatched_left, unmatched_right)."""
    pairs, used = [], set()
    for a in left:
        for j, b in enumerate(right):
            if j not in used and key(a) == key(b):
                pairs.append((a, b))
                used.add(j)
                break
    ml = {a for a, _ in pairs}
    return pairs, [a for a in left if a not in ml], [b for j, b in enumerate(right) if j not in used]


def read_adj():
    if not os.path.exists(ADJ):
        return {}
    with open(ADJ, encoding="utf-8-sig") as f:
        return {(r["pdb_id"], r["comparison"], r["side"], r["name"]): r for r in csv.DictReader(f)}


def prepare():
    data, old = load(), read_adj()
    rows = []
    for p, d in data.items():
        comparisons = [("A_vs_B", ("A", d["A"]), ("B", d["B"]))] + [
            (f"{s}_vs_B", (s, items), ("B", d["B"])) for s, items in d["systems"].items()]
        for comp, (ln, left), (rn, right) in comparisons:
            _, ul, ur = auto_match(left, right)
            for side, items, other in ((ln, ul, ur), (rn, ur, ul)):
                for name in items:
                    prev = old.get((p, comp, side, name), {})
                    rows.append({"item": d["item"], "pdb_id": p, "comparison": comp, "side": side,
                                 "name": name, "other_side_unmatched": "; ".join(other),
                                 "same_as": prev.get("same_as", ""), "consensus": prev.get("consensus", ""),
                                 "note": prev.get("note", ""), "pdbx_details": d["details"]})
    with open(ADJ, "w", newline="", encoding="utf-8-sig") as f:
        w = csv.DictWriter(f, fieldnames=list(rows[0].keys()) if rows else
                           ["item", "pdb_id", "comparison", "side", "name", "other_side_unmatched",
                            "same_as", "consensus", "note", "pdbx_details"])
        w.writeheader()
        w.writerows(rows)
    n_ab = sum(r["comparison"] == "A_vs_B" for r in rows)
    print(f"adjudication.csv: {len(rows)} unmatched items to review ({n_ab} A vs B, "
          f"{len(rows) - n_ab} extractor/LLM vs B)")


def full_match(p, comp, left_name, left, right, adj):
    """Auto matches plus adjudicated same_as pairs (decided from either side)."""
    pairs, ul, ur = auto_match(left, right)
    ur = list(ur)
    for a in list(ul):
        r = adj.get((p, comp, left_name, a))
        if r and r["same_as"].strip() in ur:
            pairs.append((a, r["same_as"].strip()))
            ur.remove(r["same_as"].strip())
            ul.remove(a)
    for b in list(ur):
        r = adj.get((p, comp, "B", b))
        if r and r["same_as"].strip() in ul:
            pairs.append((r["same_as"].strip(), b))
            ul.remove(r["same_as"].strip())
            ur.remove(b)
    return pairs, ul, ur


def ci(stat, n):
    idx = [rng.integers(0, n, n) for _ in range(N_BOOT)]
    bs = [stat(i) for i in idx]
    return [round(float(np.percentile(bs, 2.5)), 4), round(float(np.percentile(bs, 97.5)), 4)]


def prf(tp, fp, fn):
    def f(i=None):
        t, p_, n_ = (tp, fp, fn) if i is None else (tp[i], fp[i], fn[i])
        T, P, N = t.sum(), p_.sum(), n_.sum()
        prec = T / (T + P) if T + P else 0.0
        rec = T / (T + N) if T + N else 0.0
        return prec, rec, (2 * prec * rec / (prec + rec) if prec + rec else 0.0)
    tp, fp, fn = map(np.asarray, (tp, fp, fn))
    p, r, f1 = f()
    n = len(tp)
    exact = ((fp == 0) & (fn == 0)).astype(float)
    return {"tp": int(tp.sum()), "fp": int(fp.sum()), "fn": int(fn.sum()),
            "precision": round(p, 4), "precision_ci95": ci(lambda i: f(i)[0], n),
            "recall": round(r, 4), "recall_ci95": ci(lambda i: f(i)[1], n),
            "f1": round(f1, 4), "f1_ci95": ci(lambda i: f(i)[2], n),
            "entry_exact": round(float(exact.mean()), 4), "entry_exact_ci95": ci(lambda i: exact[i].mean(), n)}


def score():
    data, adj = load(), read_adj()
    missing = [k for k, r in adj.items() if k[1] == "A_vs_B" and not r["same_as"].strip()
               and r["consensus"].strip() not in ("reagent", "not_reagent")]
    ab_m, ab_a, ab_b, ab_exact = [], [], [], []
    xa = ([], [], [])
    systems = list(next(iter(data.values()))["systems"])
    vs_b = {s: ([], [], []) for s in systems}
    vs_c = {s: ([], [], []) for s in systems}
    dis = []
    for p, d in data.items():
        # ---- A vs B
        pairs, ua, ub = full_match(p, "A_vs_B", "A", d["A"], d["B"], adj)
        ab_m.append(len(pairs)); ab_a.append(len(d["A"])); ab_b.append(len(d["B"]))
        ab_exact.append(float(not ua and not ub))
        for side, items in (("A only", ua), ("B only", ub)):
            for name in items:
                r = adj.get((p, "A_vs_B", side[0], name), {})
                dis.append({"pdb_id": p, "who": side, "name": name, "consensus": r.get("consensus", ""),
                            "note": r.get("note", ""), "pdbx_details": d["details"]})
        # consensus ground truth: agreed items (B's name) + disagreements judged "reagent"
        gt_b = [b for _, b in pairs] + [b for b in ub
                                        if adj.get((p, "A_vs_B", "B", b), {}).get("consensus") == "reagent"]
        gt_a_only = [a for a in ua if adj.get((p, "A_vs_B", "A", a), {}).get("consensus") == "reagent"]
        # ---- extractor vs A (original audit)
        for k in range(3):
            xa[k].append(d["orig"][k])
        gt_n = len(gt_b) + len(gt_a_only)
        for s, items in d["systems"].items():
            # ---- system vs B
            xpairs, ux, uxb = full_match(p, f"{s}_vs_B", s, items, d["B"], adj)
            vs_b[s][0].append(len(xpairs)); vs_b[s][1].append(len(ux)); vs_b[s][2].append(len(uxb))
            # ---- system vs consensus
            x_to_b = {x: b for x, b in xpairs}
            used_x = {x for x, b in x_to_b.items() if b in gt_b}
            tp = len(used_x)
            for a in gt_a_only:  # A-only consensus items, matched to the system by name
                hit = next((x for x in items if x not in used_x and key(x) == key(a)), None)
                if hit is None:  # or via an adjudicated same_as from the system-vs-B sheet
                    hit = next((x for x in items if x not in used_x and
                                adj.get((p, f"{s}_vs_B", s, x), {}).get("same_as", "").strip() == a), None)
                if hit:
                    used_x.add(hit)
                    tp += 1
            vs_c[s][0].append(tp); vs_c[s][1].append(len(items) - len(used_x)); vs_c[s][2].append(gt_n - tp)

    ab_m, ab_a, ab_b, ab_exact = map(np.asarray, (ab_m, ab_a, ab_b, ab_exact))
    n = len(ab_m)
    psa = lambda i: 2 * ab_m[i].sum() / (ab_a[i].sum() + ab_b[i].sum())
    summary = {
        "entries": n, "unresolved_A_vs_B_disagreements": len(missing),
        "A_mentions": int(ab_a.sum()), "B_mentions": int(ab_b.sum()), "agreed_mentions": int(ab_m.sum()),
        "A_vs_B": {"positive_specific_agreement": round(float(psa(np.arange(n))), 4),
                   "psa_ci95": ci(psa, n),
                   "entry_exact_agreement": round(float(ab_exact.mean()), 4),
                   "entry_exact_ci95": ci(lambda i: ab_exact[i].mean(), n)},
        "extractor_vs_A_original": prf(*xa),
    }
    for s in systems:
        summary[f"{s}_vs_B"] = prf(*vs_b[s])
        summary[f"{s}_vs_consensus"] = (prf(*vs_c[s]) if not missing
                                        else "fill consensus for all A_vs_B rows first")
    with open(os.path.join(HERE, "agreement_summary.json"), "w") as f:
        json.dump(summary, f, indent=2)
    with open(os.path.join(HERE, "disagreements.csv"), "w", newline="", encoding="utf-8-sig") as f:
        w = csv.DictWriter(f, fieldnames=["pdb_id", "who", "name", "consensus", "note", "pdbx_details"])
        w.writeheader()
        w.writerows(dis)
    print(json.dumps(summary, indent=1))


if __name__ == "__main__":
    {"prepare": prepare, "score": score}.get(sys.argv[1] if len(sys.argv) > 1 else "", lambda: sys.exit(__doc__))()
