"""Does the structural context of the matched chain confound the
identity-vs-conditions analysis (paper Fig. `fig:retro`)?

The same protein can be crystallized alone, in a complex, as a different
oligomer, as a domain rather than the full-length chain, or with and without
a ligand. A pair of entries can therefore be 100 % identical over the
matched chain while the crystallized objects differ. This script

  1. labels the context of every entry: the polymer entities other than the
     matched chain (partners), the oligomeric state of biological assembly 1,
     and its ligands (non-polymer components minus ions, buffers,
     cryoprotectants, PEG fragments and glycans);
  2. adds to every within-target pair the alignment coverage of the two
     matched chains (aligned columns / length of the longer chain);
  3. recomputes the identity-bin curve on "same-composition" pairs: both
     entries contain a single polymer entity, have the same oligomeric
     state, and cover >= 90 % of each other (no domain-vs-full-length);
  4. splits same-UniProt pairs by context: same composition vs. different
     partners, oligomeric state, coverage or ligands.

Needs analyze_retrospective.py (pairs.csv.gz) and fetch_context.py
(entry_context.json) to have been run.

    python context_analysis.py   # -> context_summary.json, pairs_context.csv.gz
"""
import gzip
import json
import os

import numpy as np
import pandas as pd

from analyze_retrospective import ALIGNER, build_entries, clean, identity, mean_ci, spearman_ci

HERE = os.path.dirname(os.path.abspath(__file__))
MIN_COVERAGE = 90.0

# non-polymer components that say nothing about the biological ligand state:
# ions, buffers, cryoprotectants / precipitant fragments, reductants, glycans,
# unknown atoms. Everything else counts as a ligand.
ADDITIVES = {
    # ions
    "SO4", "PO4", "CL", "NA", "K", "MG", "CA", "ZN", "MN", "CU", "CU1", "FE", "FE2", "NI",
    "CO", "CD", "IOD", "BR", "NO3", "SCN", "LI", "RB", "CS", "SR", "BA", "HG", "PT", "AU",
    "YB", "SM", "GD", "TB", "PB", "F", "NH4", "AZI", "CYN", "OXY", "H2S", "NO", "CO2",
    "CMO", "PER", "O", "OH", "SE", "AL", "AF3", "ALF", "BEF", "VO4", "WO4", "MOO", "SIN",
    # buffers
    "TRS", "MES", "EPE", "IMD", "ACT", "ACY", "FMT", "CIT", "FLC", "MLI", "MLA", "TLA",
    "TAR", "BTB", "B3P", "CAC", "MPO", "PIPES", "PIG", "HEZ", "BCN", "TAM", "NHE", "CXS",
    "MOH", "OXL", "MRD", "AKG",
    # cryoprotectants, precipitant fragments, solvents, detergents
    "GOL", "EDO", "PEG", "PGE", "PG4", "1PE", "P6G", "P33", "PG0", "PG5", "2PE", "XPE",
    "12P", "15P", "PE4", "PE5", "PE8", "MPD", "IPA", "DMS", "EOH", "BU3", "PGO", "PGR",
    "DIO", "ETE", "BOG", "LDA", "C8E", "SDS", "TBU", "BME", "DTT", "DTU",
    # glycans (post-translational, not a ligand state)
    "NAG", "MAN", "BMA", "FUC", "GAL", "NDG", "SIA", "FUL", "XYP",
    # unknown
    "UNX", "UNL", "UNK",
}
# AKG, SIN, OXL are TCA/chelating anions that are almost always from the
# crystallization solution in this panel; kept as additives.


def matched_entity(target_seq, entities):
    best = None
    for ent in entities:
        if not (ent["type"] or "").startswith("polypeptide"):
            continue
        idt = identity(target_seq, clean(ent["seq"]))
        if best is None or idt > best[0]:
            best = (idt, ent)
    return best[1] if best else None


def partner_key(ent):
    """How a non-matched polymer entity is identified when comparing
    compositions: its UniProt accessions, else its polymer type."""
    if ent["uniprot"]:
        return "U:" + ",".join(sorted(ent["uniprot"]))
    return "T:" + (ent["type"] or "?")


def coverage(a, b):
    """Aligned columns / length of the longer chain, in percent."""
    if not a or not b:
        return np.nan
    c = ALIGNER.align(a, b)[0].counts()
    return 100.0 * (c.identities + c.mismatches) / max(len(a), len(b))


def build_context(E):
    with open(os.path.join(HERE, "..", "panel_10000_v2", "targets.json")) as f:
        targets = json.load(f)
    with open(os.path.join(HERE, "entry_metadata.json")) as f:
        meta = json.load(f)
    with open(os.path.join(HERE, "entry_context.json")) as f:
        ctx = json.load(f)
    rows = []
    for pid, grp in zip(E.pdb_id, E.group):
        m, c = meta[pid], ctx[pid]
        me = matched_entity(clean(targets[grp]["sequence"]), m["entities"])
        partners = tuple(sorted(partner_key(e) for e in m["entities"]
                                if me is None or e["entity"] != me["entity"]))
        ligands = frozenset(x for x in c["nonpolymer"] if x not in ADDITIVES)
        rows.append({"n_entities": len(m["entities"]), "partners": partners,
                     "olig": c["assembly1_oligomeric_state"] or "?",
                     "ligands": ligands})
    return pd.DataFrame(rows, index=E.index)


def main():
    E = build_entries()
    C = build_context(E)
    with gzip.open(os.path.join(HERE, "pairs.csv.gz"), "rt") as f:
        P = pd.read_csv(f)
    print(f"{len(P)} pairs; computing coverage ...", flush=True)
    P["coverage"] = [coverage(E.seq[i], E.seq[j]) for i, j in zip(P.i, P.j)]

    A, B = C.loc[P.i].reset_index(drop=True), C.loc[P.j].reset_index(drop=True)
    P["both_single"] = (A.n_entities.values == 1) & (B.n_entities.values == 1)
    P["same_partners"] = A.partners.values == B.partners.values
    P["alone_vs_complex"] = (A.partners.map(len).values == 0) != (B.partners.map(len).values == 0)
    P["same_olig"] = A.olig.values == B.olig.values
    P["full_coverage"] = P.coverage >= MIN_COVERAGE
    la, lb = A.ligands.values, B.ligands.values
    P["same_ligands"] = la == lb
    P["apo_vs_ligand"] = np.array([(len(x) == 0) != (len(y) == 0) for x, y in zip(la, lb)])
    # option 2 of the review: same composition, as defined in the paper
    P["same_composition"] = P.both_single & P.same_olig & P.full_coverage
    # for same-UniProt pairs, partners may be present as long as they match
    P["same_context"] = P.same_partners & P.same_olig & P.full_coverage

    with gzip.open(os.path.join(HERE, "pairs_context.csv.gz"), "wt") as f:
        P.to_csv(f, index=False)

    metrics = ["jaccard", "shared_precip", "dpH"]

    def strata(df):
        return {m: mean_ci(df[m], df.g) for m in metrics}

    labels = ["<50%", "50-70%", "70-90%", "90-98%", ">=98%"]
    P["bin"] = pd.cut(P.identity, [0, 50, 70, 90, 98, 100.01], right=False, labels=labels)
    U = P[P.uniprot_known]
    S = U[U.same_uniprot]
    Sd = S[~S.same_paper]

    def neg(df, m):
        return df[m] if m != "dpH" else -df[m]

    out = {
        "params": {"min_coverage": MIN_COVERAGE, "n_additives_excluded": len(ADDITIVES)},
        "pairs": int(len(P)),
        "fractions": {
            "all_pairs_same_composition": round(float(P.same_composition.mean()), 4),
            "ge98_same_composition": round(float(P[P.bin == ">=98%"].same_composition.mean()), 4),
            "same_uniprot_same_context": round(float(S.same_context.mean()), 4),
            "same_uniprot_alone_vs_complex": round(float(S.alone_vs_complex.mean()), 4),
            "same_uniprot_different_olig": round(float((~S.same_olig).mean()), 4),
            "same_uniprot_partial_coverage": round(float((~S.full_coverage).mean()), 4),
            "same_uniprot_apo_vs_ligand": round(float(S.apo_vs_ligand.mean()), 4),
        },
        # option 2: identity curve on same-composition pairs
        "same_composition": strata(P[P.same_composition]),
        "same_composition_by_bin": {l: strata(P[P.same_composition & (P.bin == l)]) for l in labels},
        "same_composition_bin_counts": {l: int((P.same_composition & (P.bin == l)).sum()) for l in labels},
        "different_composition_by_bin": {l: strata(P[~P.same_composition & (P.bin == l)]) for l in labels},
        "spearman_identity_vs": {
            "same_composition": {m: spearman_ci(P[P.same_composition].identity, neg(P[P.same_composition], m),
                                                P[P.same_composition].g) for m in metrics},
        },
        # option 3: same-UniProt pairs split by context
        "same_uniprot": {},
        "same_uniprot_different_paper": {},
    }
    for name, D in [("same_uniprot", S), ("same_uniprot_different_paper", Sd)]:
        sc = D.same_context
        out[name] = {
            "all": {"n": int(len(D)), **strata(D)},
            "same_context": strata(D[sc]),
            "different_context": strata(D[~sc]),
            # one factor differs, the other two held equal
            "only_partners_differ": strata(D[~D.same_partners & D.same_olig & D.full_coverage]),
            "alone_vs_complex": strata(D[D.alone_vs_complex & D.full_coverage]),
            "only_olig_differs": strata(D[D.same_partners & ~D.same_olig & D.full_coverage]),
            "only_coverage_differs": strata(D[D.same_partners & D.same_olig & ~D.full_coverage]),
            # ligand state, within otherwise identical context
            "same_context_same_ligands": strata(D[sc & D.same_ligands]),
            "same_context_apo_vs_ligand": strata(D[sc & D.apo_vs_ligand]),
            "same_context_different_ligands": strata(D[sc & ~D.same_ligands & ~D.apo_vs_ligand]),
        }

    with open(os.path.join(HERE, "context_summary.json"), "w") as f:
        json.dump(out, f, indent=2, default=lambda o: None)

    def j(d):
        return f"{d['jaccard']['mean']:.3f} {d['jaccard']['ci95']} n={d['jaccard']['n']}"
    print(json.dumps(out["fractions"], indent=1))
    print("bins same-composition:", {l: j(out["same_composition_by_bin"][l]) for l in labels})
    print("rho same-composition:", out["spearman_identity_vs"]["same_composition"])
    for name in ["same_uniprot", "same_uniprot_different_paper"]:
        print(name)
        for k, v in out[name].items():
            print(f"  {k:32s} {j(v)}  precip {v['shared_precip'].get('mean')}  dpH {v['dpH'].get('mean')}")


if __name__ == "__main__":
    main()
