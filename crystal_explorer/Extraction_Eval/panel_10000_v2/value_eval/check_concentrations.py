"""Concentration extraction on the 1,000 held-out audit entries, checked against
the source text: value, unit, and attribution to the right reagent.

Uses the extractor output frozen at audit time (../audit_1000_raw.json). Only
reagents the audit counted as correctly identified are checked (true
positives in ../annotator_A_gt/gt_1000.json).

Independent check, with no code shared with compound_extraction.py:
  1. every "number [range] unit" in the text is parsed (M, mM, uM/µM, molar,
     millimolar, %, % w/v, % v/v) and converted to M or %;
  2. the reagent's mention(s) are located by its dictionary aliases;
  3. the number-unit mention nearest to the reagent mention, in the same
     clause, is taken as the reagent's own concentration.

Classes per extracted reagent:
  correct          extracted value and unit equal the reagent's own concentration
  wrong_value      the reagent has a stated concentration, extracted value differs
  wrong_unit       same number, different unit class (M vs %)
  misattributed    extracted value equals another reagent's concentration in the text
  not_in_text      extracted value appears nowhere in the text
  missed           nothing extracted but the reagent's clause states a concentration
  none_stated      nothing extracted and nothing stated (correct)
  unlocated        the reagent mention could not be found in the text (excluded)

Everything that is not "correct" or "none_stated" goes to flagged_for_review.csv,
plus a random sample of 60 "correct" items to check the checker itself.

    python check_concentrations.py
"""
import csv
import json
import os
import random
import re
import sys

HERE = os.path.dirname(os.path.abspath(__file__))
PANEL = os.path.dirname(HERE)
CE = os.path.abspath(os.path.join(PANEL, "..", ".."))
sys.path.insert(0, CE)
from pipeline.compound_extraction import REAGENT_DICTIONARY  # noqa: E402  (aliases only, for locating mentions)

ITEM_RE = re.compile(r"'([^']+)'\s*\(([^)]*)\)")
NUM = r"(?:\d+(?:\.\d+|,\d{1,2}(?!\d))?|\.\d+)"
UNIT = (r"(?P<unit>%\s*\(?\s*(?:w\s*/\s*v|v\s*/\s*v|w\s*/\s*w|wt\s*/\s*vol|vol\s*/\s*vol)?\s*\)?|"
        r"(?:w|v)\s*/\s*v\s*%|millimolar|micromolar|molar|MM|mM|uM|UM|µM|μM|M)(?![A-Za-z])")
CONC_RE = re.compile(rf"(?<![\w.])(?P<a>{NUM})\s*(?:%\s*)?(?:(?:-|–|~|to)\s*(?P<b>{NUM}))?\s*{UNIT}")
CLAUSE_SPLIT = re.compile(r";|\n")


def to_num(s):
    return float(s.replace(",", "."))


def canon(a, b, unit):
    """(low, high, class) in M or %."""
    u = unit.lower().replace(" ", "")
    lo, hi = to_num(a), to_num(b) if b else to_num(a)
    if u.startswith("%") or u.endswith("%"):
        return lo, hi, "%"
    f = {"m": 1, "molar": 1, "mm": 1e-3, "millimolar": 1e-3, "um": 1e-6, "µm": 1e-6, "μm": 1e-6, "micromolar": 1e-6}
    k = {"M": "m", "mM": "mm", "MM": "mm", "UM": "um"}.get(unit, u)
    return lo * f[k], hi * f[k], "M"


def parse_extracted(conc):
    m = re.match(rf"\s*({NUM})(?:-({NUM}))?\s*(%|mM|M)", conc)
    if not m:
        return None
    return canon(m.group(1), m.group(2), m.group(3))


def same(x, y):
    return x[2] == y[2] and all(abs(p - q) <= 1e-6 * max(1, abs(q)) + 1e-9 for p, q in zip(x[:2], y[:2]))


def same_number(x, y):
    return abs(x[0] - y[0]) < 1e-9 or abs(x[1] - y[1]) < 1e-9


def aliases(name):
    al = {name.lower()}
    al |= {a.lower() for a in REAGENT_DICTIONARY.get(name, [])}
    return sorted(al, key=len, reverse=True)


def main():
    audit = {e["pdb_id"]: e for e in json.load(open(os.path.join(PANEL, "audit_1000_raw.json")))}
    gt = json.load(open(os.path.join(PANEL, "annotator_A_gt", "gt_1000.json")))
    rows = []
    for p, e in audit.items():
        text = (e["pdbx_details"] or "").strip().strip("'").strip(";")
        low = text.lower()
        concs = [(m.start(), m.end(), canon(m.group("a"), m.group("b"), m.group("unit")), m.group(0))
                 for m in CONC_RE.finditer(text)]
        clause_id = []
        cid = 0
        for ch_i in range(len(text) + 1):
            clause_id.append(cid)
            if ch_i < len(text) and CLAUSE_SPLIT.match(text, ch_i):
                cid += 1
        correct_ids = {x for x in gt[p]["extractor"] if x in gt[p]["reagents"]}
        seen = set()
        for name, conc in ITEM_RE.findall(e["extracted_compound_string"] or ""):
            if name in seen or name not in correct_ids:
                continue
            seen.add(name)
            spans = []
            for a in aliases(name):
                for m in re.finditer(r"(?<![a-z])" + re.escape(a) + r"(?![a-z])", low):
                    spans.append((m.start(), m.end()))
            ex = parse_extracted(conc) if conc.strip().lower() != "nan" else None
            if not spans:
                cls, own = "unlocated", None
            else:
                # nearest number-unit in the same clause as any mention of the reagent
                best = None
                for s0, s1 in spans:
                    for c0, c1, val, raw in concs:
                        if clause_id[c0] != clause_id[s0] and clause_id[c1 - 1] != clause_id[s1 - 1]:
                            continue
                        d = s0 - c1 if c1 <= s0 else (c0 - s1 if c0 >= s1 else 0)
                        if best is None or d < best[0]:
                            best = (d, val, raw)
                own = best
                if ex is None:
                    cls = "missed" if own else "none_stated"
                elif own and same(ex, own[1]):
                    cls = "correct"
                elif own and same_number(ex, own[1]) and ex[2] != own[1][2]:
                    cls = "wrong_unit"
                elif any(same(ex, v) for _, _, v, _ in concs):
                    cls = "misattributed" if own else "misattributed_no_own"
                elif own:
                    cls = "wrong_value"
                else:
                    cls = "not_in_text"
            rows.append({"pdb_id": p, "reagent": name, "extracted": conc,
                         "reagent_own_concentration_in_text": own[2] if own else "",
                         "auto_class": cls, "review_verdict": "", "note": "", "pdbx_details": text})

    with open(os.path.join(HERE, "concentration_items.csv"), "w", newline="", encoding="utf-8") as f:
        w = csv.DictWriter(f, fieldnames=list(rows[0].keys()))
        w.writeheader()
        w.writerows(rows)
    flagged = [r for r in rows if r["auto_class"] not in ("correct", "none_stated", "unlocated")]
    random.seed(7)
    sample = random.sample([r for r in rows if r["auto_class"] == "correct"], 60)
    with open(os.path.join(HERE, "flagged_for_review.csv"), "w", newline="", encoding="utf-8") as f:
        w = csv.DictWriter(f, fieldnames=["set"] + list(rows[0].keys()))
        w.writeheader()
        for r in flagged:
            w.writerow({"set": "flagged", **r})
        for r in sample:
            w.writerow({"set": "check_of_correct", **r})
    from collections import Counter
    c = Counter(r["auto_class"] for r in rows)
    print(len(rows), "correctly identified reagents checked:", dict(c.most_common()))
    print(f"flagged {len(flagged)} + 60 sampled 'correct' -> flagged_for_review.csv")


if __name__ == "__main__":
    main()
