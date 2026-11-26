"""Why each of the 167 held-out false negatives was missed, from annotator A's
notes (../annotations/). Categories:

  T  typo, abbreviation or unusual format of a reagent the dictionary has
     ("NA BR", "CdCl", "HEPS", "Sodiumcitrate", "PEG3.4K")
  R  crystallization reagent with no dictionary entry
     (succinic acid, guanidine, K3PO4, detergents, rare buffers)
  L  ligand, inhibitor, substrate or soaking compound
     (compound codes, drugs, nucleotides, sugars, metabolites)
  M  a related but different compound was reported instead
     (the FN half of an FP+FN pair: PEG MME read as PEG, triglycine as glycylglycine)

Classified by Claude from the notes; to be checked by annotator A.

    python fn_categories.py   # -> fn_categories.csv + counts
"""
import csv
import glob
import json
import os
from collections import Counter

from corrections import FIX

HERE = os.path.dirname(os.path.abspath(__file__))

M = {("5V24", "PEG MME 2000"), ("7L67", "Fuc-(alpha1-3)-GlcNAc"), ("6C4D", "triglycine"),
     ("5UBG", "PR-ATP"), ("7FCZ", "triglycine"), ("6DYT", "pyridoxal 5'-phosphate"),
     ("5GTC", "potassium cacodylate"), ("2I9T", "PEG MME 5000"), ("6ECG", "pyridoxal 5'-phosphate"),
     ("6MB1", "myristoyl-CoA"), ("3MIA", "triglycine"), ("8FY0", "tert-butanediol"),
     ("7UJR", "bis-tris propane"), ("1QI2", "2-fluoro-cellotriose"), ("4QKQ", "PEG MME 550"),
     ("3QN6", "deaza-PLP"), ("1G9M", "PEG MME 5000"), ("3IPK", "PEG MME 2000"),
     ("3P3N", "2-oxoglutarate")}
T = {"2CK2", "4TQC", "5FU0", "1YDG", "6MSC", "1DP6", "5QRP", "4BIE", "8V6Z", "1HVX", "29MX", "5TDE",
     "5CQW", "2ZVI", "9Q31", "1ECX", "9SDS", "5I78", "8V6X", "2XQA", "8P06", "1GZT", "3H5X", "4Z7H",
     "5ULW", "3HGG", "4YEJ", "1E4D", "5QRM", "4BVP", "1KEH", "7A3F", "1DWP", "4X6N", "9F8C", "5QRH",
     "3LX1", "1E3U", "9WNA", "2XQ5", "1A5A", "8C8C", "8K00", "4ZEK", "1LSX", "3IQ7", "4CIW", "5LZ2",
     "2PLF", "8V79", "4BCD", "5N47", "4B4G", "6QFT", "6I8O"}
T_ITEMS = {("3LK8", "lithium chloride"), ("3QN6", "triethanolamine"), ("6HJ3", "sodium thiocyanate"),
           ("1M9A", "DTT"), ("1M9A", "trehalose"), ("4EFS", "PEG MME 5000")}
L = {("8AD7", "D-tryptophan"), ("2Y33", "FG2"), ("2Y33", "GSNO"), ("2YPH", "Rp-cAMPS"),
     ("1R51", "8-azaxanthine"), ("6R5F", "DHP77"), ("5GK5", "docosahexaenoic acid"), ("4JC5", "CMA"),
     ("6DER", "thiamine diphosphate"), ("6DER", "metosulam"), ("5YC8", "NMS"), ("1LRM", "dihydrobiopterin"),
     ("3LK8", "AHA"), ("3P5F", "Man5"), ("5VD7", "RAC-IV-098"), ("6YCT", "p-vanillin"),
     ("6HQP", "guaiacol"), ("1RKS", "ribose"), ("6MB1", "IMP-1002"), ("6EKX", "1-naphthol"),
     ("8W0H", "HGN-1195"), ("2D3N", "maltohexaose"), ("3ZBS", "Sp-cAMPS"), ("3VV5", "AEC"),
     ("8W0L", "HGN-1196"), ("3M7W", "dehydroquinic acid"), ("6XO3", "alpha-ketoglutarate"),
     ("7PNA", "12-methoxylauric acid"), ("7TFP", "2-ketoglutarate"), ("3T9F", "5-IP7"),
     ("8E0Y", "Doxime"), ("2ZCR", "BPH-698"), ("4BVO", "POM1"), ("3UG2", "gefitinib"), ("1FXO", "TMP"),
     ("1OCQ", "cellobio-derived isofagomine"), ("1IVB", "BANA105"), ("2VS5", "UDP-galactose"),
     ("2VA3", "ddCTP"), ("6DT4", "cAMP"), ("3R60", "ascorbic acid"), ("1G93", "Hg-UDP-galactose"),
     ("5VCW", "pelitinib"), ("4N37", "methyl mannoside"), ("5T2T", "L742001"), ("9CD8", "HGN-1310"),
     ("14CQ", "phenylalanine"), ("5SCA", "phenylalanine")}


def category(p, item):
    if (p, item) in M:
        return "M"
    if (p, item) in L:
        return "L"
    if (p, item) in T_ITEMS or (p in T and FIX[p][0] == []):
        return "T"
    return "R"


def main():
    ann = {a["pdb_id"]: a for f in sorted(glob.glob(os.path.join(HERE, "..", "annotations", "*.json")))
           for a in json.load(open(f))}
    rows = [{"pdb_id": p, "missed_reagent": a, "category": category(p, a),
             "annotator_note": ann[p]["ground_truth_notes"]}
            for p, (_, adds) in FIX.items() for a in adds]
    assert len(rows) == sum(a["fn"] for a in ann.values()), (len(rows), sum(a["fn"] for a in ann.values()))
    with open(os.path.join(HERE, "fn_categories.csv"), "w", newline="", encoding="utf-8") as f:
        w = csv.DictWriter(f, fieldnames=list(rows[0].keys()))
        w.writeheader()
        w.writerows(rows)
    c = Counter(r["category"] for r in rows)
    n = len(rows)
    for k, name in (("T", "typo / abbreviation / format"), ("R", "reagent not in dictionary"),
                    ("L", "ligand / soaking compound"), ("M", "related compound reported instead")):
        print(f"{k} {name:36s} {c[k]:4d}  {100 * c[k] / n:5.1f}%")
    print(f"  total {n}")
    print("R items:", sorted({r['missed_reagent'] for r in rows if r['category'] == 'R'}))


if __name__ == "__main__":
    main()
