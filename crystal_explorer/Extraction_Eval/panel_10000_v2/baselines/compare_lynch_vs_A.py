"""Score the Lynch et al. (2020) parser against annotator A's ground truth on
the 1,000 held-out audit entries, next to our extractor's audited scores.

Inputs
  ../annotator_A_gt/gt_1000.json        A's ground truth as reagent lists (build_gt.py)
  lynch2020_all/predictions.jsonl       every name the parser returns (run_lynch2020.py)
  lynch2020_recognized/predictions.jsonl  only dictionary-recognized names
  lynch_name_matches.csv                name equivalences decided for this comparison
Output
  lynch_vs_A_summary.json, lynch_vs_A_entries.csv

Matching, applied as in A's audit of our extractor:
  - automatic: same normalized name (key() of ../second_annotator/score_agreement.py);
  - lynch_name_matches.csv, within the same entry:
      same / family        the parser's name is the same reagent or a name of the
                           right chemical family (A credited ours the same way);
      mixture_component    A counts a named mix (Tacsimate, SPG, Na/K phosphate...)
                           as one reagent and the parser returns its components:
                           one component is credited, the others are not FPs;
      fused                one string naming two or more reagents that were not
                           split: credits at most one reagent (strict) or each
                           reagent it names (lenient).
  Anything else the parser returns is a false positive, and any ground-truth reagent
  left unmatched is a false negative.

Views of the Lynch output
  parser_all_strict     every returned name; fused strings credit one reagent
  parser_all_lenient    same, fused strings credit every reagent they name
  parser_recognized     only names in their compound dictionary
  database_view         what their curated database holds: recognized names for
                        "sensible" entries, nothing for the others

Subsets: entries whose text contains "cryo", "soak" or "reservoir solution",
where the parser deliberately drops text by design, against all other entries.

95% CIs: bootstrap over entries, 5,000 resamples (as in Table I).
"""
import csv
import json
import os
import sys

import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))
PANEL = os.path.dirname(HERE)
sys.path.insert(0, os.path.join(PANEL, "second_annotator"))
from score_agreement import auto_match  # noqa: E402

N_BOOT = 5000
rng = np.random.default_rng(2034)
SCOPE_WORDS = ("cryo", "soak", "reservoir solution")


def load_jsonl(path):
    with open(path) as f:
        return {r["pdb_id"]: r for r in map(json.loads, f)}


def load_matches():
    eq = {}
    with open(os.path.join(HERE, "lynch_name_matches.csv"), encoding="utf-8") as f:
        for r in csv.DictReader(f):
            eq.setdefault(r["lynch_name"], []).append((r["gt_name"], r["kind"]))
    return eq


def score_entry(pred, gt, eq, lenient):
    """(tp, fp, fn, details) for one entry."""
    pairs, ul, ug = auto_match(pred, gt)
    tp = len(pairs)
    used = {}  # lynch name -> number of GT items it has been credited with
    absorbed = set()
    matched_g = set()
    # specific equivalences first, then mixture components, then fused strings
    for kinds in (("same", "family"), ("mixture_component",), ("fused",)):
        for g in ug:
            if g in matched_g:
                continue
            for l in ul:
                options = [k for (gname, k) in eq.get(l, []) if gname == g and k in kinds]
                if not options:
                    continue
                kind = options[0]
                if l in used and not (kind == "fused" and lenient):
                    continue
                used[l] = used.get(l, 0) + 1
                matched_g.add(g)
                tp += 1
                if kind == "mixture_component":
                    # the mix's other components are not false positives
                    for l2 in ul:
                        if l2 != l and any(gname == g and k == "mixture_component"
                                           for gname, k in eq.get(l2, [])):
                            absorbed.add(l2)
                break
    fp_items = [l for l in ul if l not in used and l not in absorbed]
    fn_items = [g for g in ug if g not in matched_g]
    return tp, len(fp_items), len(fn_items), fp_items, fn_items


def prf(tp, fp, fn):
    tp, fp, fn = map(np.asarray, (tp, fp, fn))

    def stat(i):
        T, P, N = tp[i].sum(), fp[i].sum(), fn[i].sum()
        p = T / (T + P) if T + P else 0.0
        r = T / (T + N) if T + N else 0.0
        return p, r, (2 * p * r / (p + r) if p + r else 0.0), float(((fp[i] == 0) & (fn[i] == 0)).mean())

    n = len(tp)
    point = stat(np.arange(n))
    boots = np.array([stat(rng.integers(0, n, n)) for _ in range(N_BOOT)])
    out = {"entries": n, "tp": int(tp.sum()), "fp": int(fp.sum()), "fn": int(fn.sum())}
    for k, name in enumerate(("precision", "recall", "f1", "entry_exact")):
        out[name] = round(point[k], 4)
        out[name + "_ci95"] = [round(float(np.percentile(boots[:, k], 2.5)), 4),
                               round(float(np.percentile(boots[:, k], 97.5)), 4)]
    return out


def main():
    gt = json.load(open(os.path.join(PANEL, "annotator_A_gt", "gt_1000.json")))
    audit = {e["pdb_id"]: e for e in json.load(open(os.path.join(PANEL, "audit_1000_raw.json")))}
    lall = load_jsonl(os.path.join(HERE, "lynch2020_all", "predictions.jsonl"))
    lrec = load_jsonl(os.path.join(HERE, "lynch2020_recognized", "predictions.jsonl"))
    eq = load_matches()
    ids = sorted(gt)

    views = {
        "parser_all_strict": (lambda p: lall[p]["reagents"], False),
        "parser_all_lenient": (lambda p: lall[p]["reagents"], True),
        "parser_recognized": (lambda p: lrec[p]["reagents"], False),
        "database_view": (lambda p: lrec[p]["reagents"] if lrec[p]["sensible"] else [], False),
    }
    per = {v: {} for v in views}
    for v, (get, lenient) in views.items():
        for p in ids:
            per[v][p] = score_entry(get(p), gt[p]["reagents"], eq, lenient)
    # ours: A's audited counts
    ours = {p: (gt[p]["tp"], gt[p]["fp"], gt[p]["fn"]) for p in ids}

    in_scope = {p: any(w in (audit[p]["pdbx_details"] or "").lower() for w in SCOPE_WORDS) for p in ids}
    subsets = {"all": ids,
               "text_without_cryo_soak_reservoir": [p for p in ids if not in_scope[p]],
               "text_with_cryo_soak_reservoir": [p for p in ids if in_scope[p]]}

    summary = {"ground_truth_reagents": sum(len(gt[p]["reagents"]) for p in ids),
               "lynch_sensible_entries": sum(lrec[p]["sensible"] for p in ids),
               "lynch_empty_entries": sum(not lall[p]["reagents"] for p in ids)}
    for sname, sids in subsets.items():
        block = {"entries": len(sids),
                 "ours_extractor": prf(*zip(*[ours[p] for p in sids]))}
        for v in views:
            block["lynch_" + v] = prf(*zip(*[per[v][p][:3] for p in sids]))
        summary[sname] = block
    # their database only describes sensible entries: score those on their own too
    sens = [p for p in ids if lrec[p]["sensible"]]
    summary["sensible_entries_only"] = {"entries": len(sens),
                                        "ours_extractor": prf(*zip(*[ours[p] for p in sens])),
                                        "lynch_database_view": prf(*zip(*[per["database_view"][p][:3] for p in sens]))}

    with open(os.path.join(HERE, "lynch_vs_A_summary.json"), "w") as f:
        json.dump(summary, f, indent=2)
    with open(os.path.join(HERE, "lynch_vs_A_entries.csv"), "w", newline="", encoding="utf-8") as f:
        w = csv.writer(f)
        w.writerow(["pdb_id", "scope_words_in_text", "lynch_sensible", "gt", "lynch_all",
                    "strict_tp", "strict_fp", "strict_fn", "strict_fp_items", "strict_fn_items",
                    "ours_tp", "ours_fp", "ours_fn", "pdbx_details"])
        for p in ids:
            tp, fp, fn, fpi, fni = per["parser_all_strict"][p]
            w.writerow([p, in_scope[p], lrec[p]["sensible"], "; ".join(gt[p]["reagents"]),
                        "; ".join(lall[p]["reagents"]), tp, fp, fn, "; ".join(fpi), "; ".join(fni),
                        *ours[p], audit[p]["pdbx_details"]])

    def row(d):
        return (f"P {d['precision']:.3f} [{d['precision_ci95'][0]:.3f}-{d['precision_ci95'][1]:.3f}]  "
                f"R {d['recall']:.3f} [{d['recall_ci95'][0]:.3f}-{d['recall_ci95'][1]:.3f}]  "
                f"F1 {d['f1']:.3f}  exact {d['entry_exact']:.3f}  (tp {d['tp']}, fp {d['fp']}, fn {d['fn']})")
    for sname in subsets:
        print(f"\n== {sname} ({summary[sname]['entries']} entries)")
        for k, d in summary[sname].items():
            if isinstance(d, dict):
                print(f"  {k:28s} {row(d)}")
    print(f"\n== sensible entries only ({len(sens)})")
    for k in ("ours_extractor", "lynch_database_view"):
        print(f"  {k:28s} {row(summary['sensible_entries_only'][k])}")


if __name__ == "__main__":
    main()
