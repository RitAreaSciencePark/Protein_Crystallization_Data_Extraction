"""Build the blinded second-annotator package for the held-out audit.

Draws 200 of the 1,000 audited entries (seed 2032) and writes:

  for_annotator_B/annotation_sheet.csv   PDB ID + raw pdbx_details only, with an
                                         empty `reagents` column for annotator B.
                                         No extractor output, no annotator-A notes.
  for_annotator_B/GUIDELINES.md          copied from ./GUIDELINES.md
  annotator_A_structured.csv             annotator A's original ground truth for
                                         the same 200 entries as a reagent list,
                                         to be CHECKED by annotator A (see below)
  sample_meta.json                       seed, ids, SHA-256 of inputs

Annotator A's original annotation is free text plus TP/FP/FN counts. For an
entry with no errors (fp = fn = 0) A's ground truth is exactly the
extractor's list, so the structured list is filled automatically (and
checked against the TP count). For entries with errors, the draft is the
extractor's list plus every quoted missed name from the notes, flagged
needs_review = yes: A must remove the false positives and fix the names by
hand before scoring. Only this file should be edited by A; it must be
finalized BEFORE A sees annotator B's sheet.

    python build_package.py
"""
import csv
import glob
import hashlib
import json
import os
import random
import re
import shutil

HERE = os.path.dirname(os.path.abspath(__file__))
PANEL = os.path.dirname(HERE)
SEED = 2032
N = 200
ITEM_RE = re.compile(r"'([^']+)'\s*\(")


def sha(path):
    with open(path, "rb") as f:
        return hashlib.sha256(f.read()).hexdigest()


def main():
    audit_path = os.path.join(PANEL, "audit_1000_raw.json")
    with open(audit_path) as f:
        audit = json.load(f)
    ann = {}
    for p in sorted(glob.glob(os.path.join(PANEL, "annotations", "batch_*_annotations.json"))):
        for a in json.load(open(p)):
            ann[a["pdb_id"]] = a
    assert len(audit) == 1000 and len(ann) == 1000, (len(audit), len(ann))

    sample = random.Random(SEED).sample(audit, N)
    out_b = os.path.join(HERE, "for_annotator_B")
    os.makedirs(out_b, exist_ok=True)

    # ---- annotator B: blinded sheet -------------------------------------
    with open(os.path.join(out_b, "annotation_sheet.csv"), "w", newline="", encoding="utf-8-sig") as f:
        w = csv.writer(f)
        w.writerow(["item", "pdb_id", "pdbx_details", "reagents", "comments"])
        for k, e in enumerate(sample, 1):
            w.writerow([k, e["pdb_id"], e["pdbx_details"].strip("'").strip(";").strip(), "", ""])
    shutil.copy(os.path.join(HERE, "GUIDELINES.md"), os.path.join(out_b, "GUIDELINES.md"))

    # ---- annotator A: structured draft ----------------------------------
    n_review = 0
    with open(os.path.join(HERE, "annotator_A_structured.csv"), "w", newline="", encoding="utf-8-sig") as f:
        w = csv.writer(f)
        w.writerow(["item", "pdb_id", "A_reagents", "needs_review", "tp", "fp", "fn",
                    "extractor_reagents", "A_original_notes", "pdbx_details"])
        for k, e in enumerate(sample, 1):
            a = ann[e["pdb_id"]]
            ext = list(dict.fromkeys(ITEM_RE.findall(e["extracted_compound_string"] or "")))
            draft = list(ext)
            review = bool(a["fp"] or a["fn"])
            if review:
                # quoted names in the notes are, by A's convention, the missed source strings
                for q in re.findall(r"'([^']+)'", a["ground_truth_notes"]):
                    if q not in draft:
                        draft.append(f"{q} [FN?]")
            elif len(ext) != a["tp"]:
                review = True  # extractor list and A's TP count disagree: check by hand
            n_review += review
            w.writerow([k, e["pdb_id"], "; ".join(draft), "yes" if review else "no",
                        a["tp"], a["fp"], a["fn"], "; ".join(ext), a["ground_truth_notes"],
                        e["pdbx_details"]])

    meta = {"seed": SEED, "n": N, "pdb_ids": [e["pdb_id"] for e in sample],
            "audit_1000_raw_sha256": sha(audit_path),
            "extractor_sha256_at_audit": json.load(open(os.path.join(PANEL, "panel_meta.json")))[
                "compound_extraction_sha256"],
            "entries_needing_A_review": n_review}
    with open(os.path.join(HERE, "sample_meta.json"), "w") as f:
        json.dump(meta, f, indent=2)
    print(f"{N} entries sampled; {n_review} need annotator-A review in annotator_A_structured.csv")


if __name__ == "__main__":
    main()
