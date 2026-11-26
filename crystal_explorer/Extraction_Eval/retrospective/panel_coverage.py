"""Diversity of the held-out 10,000-entry panel: how many unique proteins,
organisms, depositors/publications and condition types it covers.

    python panel_coverage.py   # needs entry_metadata.json + entry_metadata_extra.json
                               # -> panel_coverage_stats.json
"""
import json
import os
import re
from collections import Counter

import numpy as np

from analyze_retrospective import PANEL_DIR, clean, identity, parse_reagents

HERE = os.path.dirname(os.path.abspath(__file__))


def method_class(m):
    m = (m or "").lower()
    if not m or m in ("?", "."):
        return "not reported"
    if "hanging" in m:
        return "vapor diffusion, hanging drop"
    if "sitting" in m:
        return "vapor diffusion, sitting drop"
    if "vapor" in m or "vapour" in m:
        return "vapor diffusion, unspecified"
    if "batch" in m or "microbatch" in m:
        return "batch / microbatch"
    if "lipidic" in m or "lcp" in m or "cubic" in m:
        return "lipidic cubic phase"
    if "dialysis" in m:
        return "dialysis"
    if "evaporation" in m:
        return "evaporation"
    if "counter" in m or "diffusion" in m:
        return "other diffusion"
    return "other"


def main():
    panel = json.load(open(os.path.join(PANEL_DIR, "panel_10000_raw.json")))
    targets = json.load(open(os.path.join(PANEL_DIR, "targets.json")))
    meta = json.load(open(os.path.join(HERE, "entry_metadata.json")))
    extra = json.load(open(os.path.join(HERE, "entry_metadata_extra.json")))

    uniprots, seqs, organisms, kingdoms = set(), set(), Counter(), Counter()
    authors, last_authors, first_authors, pubs = set(), set(), set(), set()
    reag_count, reag_sets, precip_count, methods = Counter(), set(), Counter(), Counter()
    pHs, temps, years, lengths = [], [], [], []
    no_pub = 0
    for e in panel:
        pid = e["pdb_id"].upper()
        m, x = meta[pid], extra[pid]
        tseq = clean(targets[e["protein"]]["sequence"])
        best = max(((identity(tseq, clean(ent["seq"])), ent) for ent in m["entities"]
                    if (ent["type"] or "").startswith("polypeptide")), key=lambda t: t[0])[1]
        uniprots.update(best["uniprot"])
        seqs.add(best["seq"])
        lengths.append(len(best["seq"]))
        for o in x["organisms"].get(best["entity"], []) or ["unknown"]:
            organisms[o] += 1
        a = [n.strip() for n in x["authors"] if n and n.strip()]
        authors.update(a)
        if a:
            first_authors.add(a[0])
            last_authors.add(a[-1])
        if m["pubmed"]:
            pubs.add(m["pubmed"])
        else:
            no_pub += 1
        reag, prec = parse_reagents(e["extracted_compound_string"])
        reag_count.update(reag)
        if reag:
            reag_sets.add(reag)
        precip_count.update(prec)
        methods[method_class((x["method"] or [None])[0])] += 1
        if m["pH"] is not None:
            pHs.append(m["pH"])
        if m["temp"] is not None:
            temps.append(m["temp"])
        if m["deposit_date"]:
            years.append(int(m["deposit_date"][:4]))

    def precip_family(p):
        if p.startswith("PEG"):
            return p
        if any(w in p for w in ("mpd", "isopropanol", "ethanol", "hexanediol", "jeffamine",
                                "pentaerythritol", "butanol", "methanol", "butanediol")):
            return "organic (MPD, alcohols, ...)"
        if p == "tacsimate":
            return "salt (>= 0.8 M)"
        return "salt (>= 0.8 M)"

    fam = Counter()
    for p, n in precip_count.items():
        fam[precip_family(p)] += n

    pHs, temps = np.array(pHs, float), np.array(temps, float)
    org_top = organisms.most_common(8)
    stats = {
        "entries": len(panel),
        "targets_in_panel": len({e["protein"] for e in panel}),
        "targets_run_ok": 683,
        "clusters_30pct_eligible": 948,
        "clusters_30pct_total": 25325,
        "proteins": {
            "unique_uniprot_accessions": len(uniprots),
            "unique_matched_chain_sequences": len(seqs),
            "unique_source_organisms": len(organisms),
            "top_organisms": org_top,
            "chain_length_median": int(np.median(lengths)),
            "chain_length_range": [int(min(lengths)), int(max(lengths))],
        },
        "depositors": {
            "unique_author_names": len(authors),
            "unique_last_authors (proxy for labs)": len(last_authors),
            "unique_first_authors": len(first_authors),
            "unique_primary_publications (PubMed)": len(pubs),
            "entries_without_pubmed": no_pub,
            "deposition_years": [min(years), max(years)],
            "deposition_year_median": int(np.median(years)),
        },
        "conditions": {
            "distinct_reagents": len(reag_count),
            "distinct_reagent_combinations": len(reag_sets),
            "distinct_precipitants": len(precip_count),
            "precipitant_families (entry mentions)": dict(fam.most_common()),
            "pH": {"n": int(len(pHs)), "min": float(pHs.min()), "max": float(pHs.max()),
                   "p5_p50_p95": [float(np.percentile(pHs, q)) for q in (5, 50, 95)],
                   "distinct_values": int(len(np.unique(pHs)))},
            "temperature_K": {"n": int(len(temps)), "p5_p50_p95": [float(np.percentile(temps, q)) for q in (5, 50, 95)],
                              "share_277K": round(float(np.mean(np.abs(temps - 277) <= 1.5)), 3),
                              "share_291_298K": round(float(np.mean((temps >= 290.5) & (temps <= 298.5))), 3)},
            "crystallization_method": dict(methods.most_common()),
        },
    }
    with open(os.path.join(HERE, "panel_coverage_stats.json"), "w") as f:
        json.dump(stats, f, indent=2)
    print(json.dumps(stats, indent=1))


if __name__ == "__main__":
    main()
