"""Annotator A's ground truth for all 1,000 held-out audit entries as reagent
lists (gt_1000.json), reconstructed from the original audit:

- error-free entries (fp = fn = 0): the extractor's own list, which A accepted;
- entries with errors: the extractor's list edited per corrections.py, which
  was written from A's notes.

Every entry is checked against A's audit counts: the list must hold tp + fn
items, and after the edits the extractor's list must contain exactly tp of
them and fp others. Rerunning the original scoring on this file therefore
reproduces the audited precision and recall, apart from the documented
exceptions in corrections.COUNT_EXCEPTIONS.

    python build_gt.py
"""
import glob
import json
import os
import re

from corrections import COUNT_EXCEPTIONS, FIX

HERE = os.path.dirname(os.path.abspath(__file__))
PANEL = os.path.dirname(HERE)
ITEM_RE = re.compile(r"'([^']+)'\s*\(")


def main():
    ann = {a["pdb_id"]: a for p in sorted(glob.glob(os.path.join(PANEL, "annotations", "batch_*_annotations.json")))
           for a in json.load(open(p))}
    audit = json.load(open(os.path.join(PANEL, "audit_1000_raw.json")))
    gt, problems = {}, []
    for e in audit:
        p, a = e["pdb_id"], ann[e["pdb_id"]]
        ext = list(dict.fromkeys(ITEM_RE.findall(e["extracted_compound_string"] or "")))
        if a["fp"] or a["fn"] or p in FIX:
            if p not in FIX:
                problems.append(f"{p}: has errors but no entry in corrections.py")
                continue
            remove, add = FIX[p]
            missing = [r for r in remove if r not in ext]
            if missing:
                problems.append(f"{p}: items to remove not in extractor list: {missing}")
            items = [x for x in ext if x not in remove] + add
            tp_implied, fp_implied, fn_implied = len(ext) - len(remove), len(remove), len(add)
        else:
            items = ext
            tp_implied, fp_implied, fn_implied = len(ext), 0, 0
        if (tp_implied, fp_implied, fn_implied) != (a["tp"], a["fp"], a["fn"]) and p not in COUNT_EXCEPTIONS:
            problems.append(f"{p}: implied tp/fp/fn {tp_implied}/{fp_implied}/{fn_implied} "
                            f"!= audited {a['tp']}/{a['fp']}/{a['fn']}")
        gt[p] = {"reagents": items, "extractor": ext, "tp": a["tp"], "fp": a["fp"], "fn": a["fn"],
                 "source": "corrections.py" if p in FIX else "extractor list accepted by A"}
    if problems:
        print("\n".join(problems))
        raise SystemExit(f"{len(problems)} problem(s); nothing written")
    with open(os.path.join(HERE, "gt_1000.json"), "w") as f:
        json.dump(gt, f, indent=1)
    n = sum(len(v["reagents"]) for v in gt.values())
    print(f"gt_1000.json: {len(gt)} entries, {n} ground-truth reagents "
          f"({sum(1 for v in gt.values() if v['source'] == 'corrections.py')} from corrections)")


if __name__ == "__main__":
    main()
