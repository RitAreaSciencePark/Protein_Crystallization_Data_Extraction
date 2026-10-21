"""Aggregate the finished audit batches (annotations/batch_NN_annotations.json)
into precision / recall / F1 / entry-level exact-match accuracy, using the
same entity-level TP/FP/FN scoring as the earlier panels. Works on a
partial audit too: it reports how many of the 20 batches are in.

Writes audit_summary.json next to this script."""
import glob
import json
import os

HERE = os.path.dirname(os.path.abspath(__file__))


def main():
    files = sorted(glob.glob(os.path.join(HERE, "annotations", "batch_*_annotations.json")))
    entries = [e for f in files for e in json.load(open(f))]
    tp = sum(e["tp"] for e in entries)
    fp = sum(e["fp"] for e in entries)
    fn = sum(e["fn"] for e in entries)
    exact = sum(1 for e in entries if e["fp"] == 0 and e["fn"] == 0)
    p = tp / (tp + fp) if tp + fp else 0.0
    r = tp / (tp + fn) if tp + fn else 0.0
    f1 = 2 * p * r / (p + r) if p + r else 0.0
    summary = {
        "batches_done": len(files), "batches_total": 20, "n_entries": len(entries),
        "total_ground_truth_mentions": tp + fn, "total_extracted_mentions": tp + fp,
        "true_positives": tp, "false_positives": fp, "false_negatives": fn,
        "precision": round(p, 4), "recall": round(r, 4), "f1": round(f1, 4),
        "entry_level_exact_match_accuracy": round(exact / len(entries), 4) if entries else 0.0,
        "entry_level_exact_match_count": exact,
    }
    with open(os.path.join(HERE, "audit_summary.json"), "w") as f:
        json.dump(summary, f, indent=2)
    print(json.dumps(summary, indent=2))


if __name__ == "__main__":
    main()
