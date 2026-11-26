"""pH and temperature text fallbacks of pipeline/pdb_sequence_search.py.

The pipeline reads _exptl_crystal_grow.pH / .temp and parses pdbx_details only
when the field is empty. Here the pipeline's own functions are called with the
structured field blanked, on all 10,000 held-out panel entries, and compared
with the deposited value (../../retrospective/entry_metadata.json).

    python check_ph_temp.py   # -> ph_temp_summary.json
"""
import json
import os
import re
import sys

HERE = os.path.dirname(os.path.abspath(__file__))
PANEL = os.path.dirname(HERE)
CE = os.path.abspath(os.path.join(PANEL, "..", ".."))
sys.path.insert(0, CE)
from pipeline.pdb_sequence_search import (  # noqa: E402
    get_ph_from_mmcif_or_details, get_temperature_from_mmcif_or_details)

AUTO_TEMP = re.compile(r"temperature\s+[0-9.]+\s*K\b", re.I)


class TextOnly:
    def __init__(self, t):
        self.t = t

    def find_value(self, tag):
        return self.t if tag == "_exptl_crystal_grow.pdbx_details" else None


def main():
    panel = json.load(open(os.path.join(PANEL, "panel_10000_raw.json")))
    meta = json.load(open(os.path.join(CE, "Extraction_Eval", "retrospective", "entry_metadata.json")))
    R = []
    for e in panel:
        t = (e["pdbx_details"] or "").strip().strip("'").strip(";")
        m = meta[e["pdb_id"].upper()]
        b = TextOnly(t)
        R.append({"ph_ref": m["pH"], "ph_fb": get_ph_from_mmcif_or_details(b),
                  "t_ref": m["temp"], "t_fb": get_temperature_from_mmcif_or_details(b),
                  "auto": bool(AUTO_TEMP.search(t)), "n_ph": len(re.findall(r"\bpH\s*[:=]?\s*\d", t, re.I))})
    both_ph = [r for r in R if r["ph_ref"] is not None and r["ph_fb"] is not None]
    agree = lambda rs, k1, k2, tol: sum(abs(r[k1] - r[k2]) <= tol for r in rs)
    s = {"entries": len(R)}
    s["pH"] = {
        "field_missing_so_fallback_used": sum(r["ph_ref"] is None for r in R),
        "fallback_finds_value_when_used": sum(r["ph_ref"] is None and r["ph_fb"] is not None for r in R),
        "implausible_when_used": sum(r["ph_ref"] is None and r["ph_fb"] is not None and not 2 <= r["ph_fb"] <= 12 for r in R),
        "validation_entries_with_deposited_pH": sum(r["ph_ref"] is not None for r in R),
        "validation_fallback_found": len(both_ph),
        "agree_exact": agree(both_ph, "ph_ref", "ph_fb", 0.001),
        "agree_within_0.1": agree(both_ph, "ph_ref", "ph_fb", 0.1001),
        "agree_within_0.5": agree(both_ph, "ph_ref", "ph_fb", 0.5001),
        "disagreements_with_several_pH_in_text": sum(1 for r in both_ph if abs(r["ph_ref"] - r["ph_fb"]) > 0.1001 and r["n_ph"] > 1),
    }
    used_t = [r for r in R if r["t_ref"] is None and r["t_fb"] is not None]
    both_t = [r for r in R if r["t_ref"] is not None and r["t_fb"] is not None and not r["auto"]]
    s["temperature"] = {
        "field_missing_so_fallback_used": sum(r["t_ref"] is None for r in R),
        "fallback_finds_value_when_used": len(used_t),
        "implausible_when_used_(outside_250-320K)": sum(not 250 <= r["t_fb"] <= 320 for r in used_t),
        "validation_texts_without_auto_appended_temperature": len(both_t),
        "agree_within_1K": agree(both_t, "t_ref", "t_fb", 1.001),
    }
    for k in ("agree_exact", "agree_within_0.1", "agree_within_0.5"):
        s["pH"][k.replace("agree", "rate")] = round(s["pH"][k] / len(both_ph), 4)
    s["pH"]["fallback_coverage"] = round(len(both_ph) / s["pH"]["validation_entries_with_deposited_pH"], 4)
    json.dump(s, open(os.path.join(HERE, "ph_temp_summary.json"), "w"), indent=2)
    print(json.dumps(s, indent=1))


if __name__ == "__main__":
    main()
