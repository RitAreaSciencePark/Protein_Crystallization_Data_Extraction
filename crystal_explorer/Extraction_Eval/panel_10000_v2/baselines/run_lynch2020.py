"""Baseline: the published parser of Lynch, Dudek & Bowman (2020, Patterns 1:100024),
run on exactly the same pdbx_details texts as our held-out audit.

We run their code, not their released database. The database is a December
2019 snapshot, and many audited entries were deposited later; running the
parser on the same text also removes any difference in input. The parser is
the upstream repository at a pinned commit, unmodified:

    git clone https://github.com/Hauptman-Woodward/crystallizationDatabase lynch
    git -C lynch checkout 42be6d1171ba1c2487e8cc4d22d86c7efd3bf84c
    python -m venv lynch_venv && lynch_venv/bin/pip install nltk
    NLTK_DATA=<dir> lynch_venv/bin/python -c "import nltk; nltk.download('punkt'); nltk.download('punkt_tab')"

    NLTK_DATA=<dir> lynch_venv/bin/python run_lynch2020.py --lynch-dir lynch --entries all

Their pipeline for one entry, reproduced exactly:
  Structure.parseDetails()   (errors caught and the entry kept, as in their parseAllDetails)
  Structure.standardizeNames()  after updateMiscDictionaries() (as in their standardizeAllNames)
Their dictionary-maintenance functions also rewrite their JSON input files;
that file writing is disabled here so the pinned clone stays untouched. The
in-memory dictionaries are updated exactly as in their code.

Two systems are written, both scored by ../second_annotator/score_agreement.py:
  lynch2020_all          every compound name the parser returns
  lynch2020_recognized   only names in their compound dictionary, i.e. what their
                         database's 0/1 compound matrix contains
plus, per entry, their "sensible" flag (every extracted name recognized and at
least one compound): only sensible entries enter their curated database.

Design differences worth reporting (they are choices, not bugs): the parser
cuts the text at "cryo" and at "soaked", and keeps only the "reservoir
solution" part when that phrase is present. It therefore does not report
cryoprotectants, soaked ligands or protein-buffer components, which our
scoring rules count.
"""
import argparse
import json
import os
import subprocess
import sys

HERE = os.path.dirname(os.path.abspath(__file__))
PANEL = os.path.dirname(HERE)
PINNED = "42be6d1171ba1c2487e8cc4d22d86c7efd3bf84c"


def clean(details):
    return (details or "").strip().strip("'").strip(";").strip()


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--lynch-dir", required=True, help="clone of Hauptman-Woodward/crystallizationDatabase")
    ap.add_argument("--entries", choices=["sample200", "all"], default="all")
    args = ap.parse_args()

    lynch = os.path.abspath(args.lynch_dir)
    commit = subprocess.run(["git", "-C", lynch, "rev-parse", "HEAD"], capture_output=True, text=True).stdout.strip()
    if commit != PINNED:
        sys.exit(f"{lynch} is at {commit}, expected the pinned commit {PINNED}")
    dirty = subprocess.run(["git", "-C", lynch, "status", "--porcelain", "--", "src", "Input"],
                           capture_output=True, text=True).stdout.strip()
    if dirty:
        sys.exit(f"the Lynch clone has local changes, refusing to run:\n{dirty}")

    with open(os.path.join(PANEL, "audit_1000_raw.json")) as f:
        audit = json.load(f)
    if args.entries == "sample200":
        with open(os.path.join(PANEL, "second_annotator", "sample_meta.json")) as f:
            ids = set(json.load(f)["pdb_ids"])
        audit = [e for e in audit if e["pdb_id"] in ids]

    # Their module resolves Input/ relative to the working directory
    os.chdir(lynch)
    sys.path.insert(0, os.path.join(lynch, "src"))
    import pdb_crystal_database as pcd  # noqa: E402
    import nltk  # noqa: E402

    pcd.writeJson = lambda *a, **k: None  # keep their input files untouched
    pcd.importNLTK()
    pcd.updateMiscDictionaries()  # in-memory part of their standardizeAllNames
    known = set(pcd.compoundDictionary.values())

    rows, n_parse_err, n_std_err = [], 0, 0
    for e in audit:
        s = pcd.Structure(e["pdb_id"], None, clean(e["pdbx_details"]), None, None, None, None, [], None)
        try:
            s.parseDetails()
        except Exception:
            n_parse_err += 1
        try:
            s.standardizeNames()
        except Exception:
            n_std_err += 1
        names = []
        for c in s.compounds[::2]:
            if c is not None and str(c).strip() and c not in names:
                names.append(c)
        recognized = [c for c in names if c in known]
        sensible = bool(names) and len(recognized) == len(names)
        rows.append({"pdb_id": e["pdb_id"], "all": names, "recognized": recognized, "sensible": sensible,
                     "raw_compounds": [str(x) if x is not None else None for x in s.compounds]})

    out_meta = {"lynch_commit": commit, "nltk_version": nltk.__version__, "entries": args.entries,
                "n": len(rows), "parse_errors": n_parse_err, "standardize_errors": n_std_err,
                "sensible_fraction": round(sum(r["sensible"] for r in rows) / len(rows), 4),
                "empty_fraction": round(sum(not r["all"] for r in rows) / len(rows), 4),
                "mean_names_all": round(sum(len(r["all"]) for r in rows) / len(rows), 3),
                "mean_names_recognized": round(sum(len(r["recognized"]) for r in rows) / len(rows), 3)}
    for system, field in (("lynch2020_all", "all"), ("lynch2020_recognized", "recognized")):
        d = os.path.join(HERE, system)
        os.makedirs(d, exist_ok=True)
        with open(os.path.join(d, "predictions.jsonl"), "w") as f:
            for r in rows:
                f.write(json.dumps({"pdb_id": r["pdb_id"], "reagents": r[field], "sensible": r["sensible"],
                                    "raw_compounds": r["raw_compounds"]}) + "\n")
        with open(os.path.join(d, "run_meta.json"), "w") as f:
            json.dump({**out_meta, "system": system, "reagents_field": field}, f, indent=2)
    print(json.dumps(out_meta, indent=2))


if __name__ == "__main__":
    main()
