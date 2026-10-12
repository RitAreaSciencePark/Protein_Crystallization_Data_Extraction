"""Compute precision/recall/F1/accuracy of compound_extraction.py against
the manually curated ground truth in ground_truth_annotations.json (60
real, randomly sampled pdbx_details entries -- see build_ground_truth_sample
in the paper's Evaluation section, Section VII). Writes eval_summary.json
next to this script."""
import json
import os
from collections import Counter

HERE = os.path.dirname(os.path.abspath(__file__))

with open(os.path.join(HERE, "ground_truth_annotations.json")) as f:
    data = json.load(f)

n = len(data)
tp = sum(d["tp"] for d in data)
fp = sum(d["fp"] for d in data)
fn = sum(d["fn"] for d in data)
conc_miss = sum(d.get("conc_miss", 0) for d in data)

precision = tp / (tp + fp)
recall = tp / (tp + fn)
f1 = 2 * precision * recall / (precision + recall)
exact_match = sum(1 for d in data if d["fp"] == 0 and d["fn"] == 0)

summary = {
    "n_entries": n,
    "total_ground_truth_mentions": tp + fn,
    "total_extracted_mentions": tp + fp,
    "true_positives": tp,
    "false_positives": fp,
    "false_negatives": fn,
    "precision": round(precision, 4),
    "recall": round(recall, 4),
    "f1": round(f1, 4),
    "entry_level_exact_match_accuracy": round(exact_match / n, 4),
    "entry_level_exact_match_count": exact_match,
    "concentration_attachment_misses": conc_miss,
}

print(json.dumps(summary, indent=2))
with open(os.path.join(HERE, "eval_summary.json"), "w") as f:
    json.dump(summary, f, indent=2)
