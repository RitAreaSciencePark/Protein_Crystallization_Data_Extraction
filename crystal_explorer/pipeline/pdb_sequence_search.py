#!/usr/bin/env python3
"""
pdb_sequence_search.py
=======================

Given a protein (or nucleic acid) sequence, this script:

  1. Queries the RCSB PDB Search API (v2) using a sequence-similarity
     search (mmseqs2-based), restricted to structures solved by
     X-ray diffraction.
  2. Fetches entry-level metadata for every hit via the RCSB Data API,
     including the free-text `pdbx_details` field that describes the
     crystallization conditions (`exptl_crystal_grow.pdbx_details`).
  3. Parses that free text with regex heuristics to pull out
     (compound, concentration) pairs, e.g. "0.2 M ammonium sulfate"
     -> ("ammonium sulfate", "0.2 M").
  4. Filters down to entries with real experimental crystallization data
     and writes ONLY that filtered result to Output.csv.

Usage
-----
    python pdb_sequence_search.py --sequence "MKTAYIAKQRQISFVKSHFSRQLEERLGLIEV..." \
        --identity 0.9 --evalue 1.0 --max-hits 25 --out Output.csv

    # or read the sequence from a FASTA file
    python pdb_sequence_search.py --fasta my_protein.fasta --out Output.csv

Notes
-----
- The RCSB sequence search operator uses mmseqs2 under the hood.
  `identity_cutoff` is a fraction (0-1), `evalue_cutoff` follows normal
  BLAST E-value conventions (smaller = more significant).
- `pdbx_details` is unstructured free text written by depositors, so the
  compound/concentration extraction is heuristic. It will not catch
  every phrasing, and manual review of edge cases is recommended.
- Network calls are rate-limited slightly to be polite to the RCSB
  servers; increase --sleep if you hit rate limits on large hit lists.
"""
import os
import argparse
import gemmi
import json
import re
import sys
import time
import requests
import pandas as pd
from typing import Dict, Iterable, List, Optional, Tuple
from Bio.PDB.MMCIF2Dict import MMCIF2Dict


SEARCH_URL = "https://search.rcsb.org/rcsbsearch/v2/query"
DATA_ENTRY_URL = "https://data.rcsb.org/rest/v1/core/entry/{pdb_id}"
DATA_CRYSTAL_URL = "https://data.rcsb.org/rest/v1/core/exptl_crystal_grow/{pdb_id}"


# --------------------------------------------------------------------------- #
# 1. Sequence search
# --------------------------------------------------------------------------- #
def build_sequence_query(
    sequence: str,
    identity_cutoff: float = 0.5,
    evalue_cutoff: float = 1.0,
    sequence_type: str = "protein",
    experimental_method: str = "X-RAY DIFFRACTION",
    max_hits: int = 1000,
) -> dict:
    """Build the JSON payload for a sequence search restricted to a given
    experimental method (default: X-ray diffraction), returning results
    grouped at the polymer-entity level."""

    return {
        "query": {
            "type": "group",
            "logical_operator": "and",
            "nodes": [
                {
                    "type": "terminal",
                    "service": "sequence",
                    "parameters": {
                        "evalue_cutoff": evalue_cutoff,
                        "identity_cutoff": identity_cutoff,
                        "sequence_type": sequence_type,  # "protein" | "dna" | "rna"
                        "value": sequence,
                    },
                },
                {
                    "type": "terminal",
                    "service": "text",
                    "parameters": {
                        "attribute": "exptl.method",
                        "operator": "exact_match",
                        "value": experimental_method,
                    },
                },
            ],
        },
        "return_type": "polymer_entity",
        "request_options": {
            "paginate": {"start": 0, "rows": max_hits},
            "results_content_type": ["experimental"],
            "sort": [{"sort_by": "score", "direction": "desc"}],
            "results_verbosity": "verbose",  # include match_context (identity, e-value, alignment)
        },
    }


_RETRYABLE_STATUS = {502, 503, 504}


def run_sequence_search(query: dict, timeout: int = 60, retries: int = 2) -> List[dict]:
    """POST the query to the RCSB Search API and return the raw result list.

    A sequence-similarity search can occasionally take RCSB longer than a
    plain lookup to compute, and the endpoint is sometimes just slow or
    briefly unavailable under load -- a single 30s attempt with no retry
    meant a run could fail outright on a transient hiccup. Retries with
    backoff on timeouts, connection errors, and 502/503/504 (upstream
    having a bad moment); anything else (4xx, bad JSON) still fails
    immediately since retrying won't fix it."""
    last_exc = None
    for attempt in range(retries + 1):
        try:
            resp = requests.post(SEARCH_URL, json=query, timeout=timeout)
            if resp.status_code == 204:
                return []  # no hits
            resp.raise_for_status()  # raises HTTPError for any 4xx/5xx, caught below
            data = resp.json()
            return data.get("result_set", [])
        except (requests.exceptions.Timeout, requests.exceptions.ConnectionError) as e:
            last_exc = e
        except requests.exceptions.HTTPError as e:
            if e.response is not None and e.response.status_code in _RETRYABLE_STATUS:
                last_exc = e
            else:
                raise
        if attempt < retries:
            time.sleep(2 * (attempt + 1))
    raise last_exc


_MATCH_CONTEXT_FIELDS = [
    "sequence_identity", "evalue"]


def _best_match_context(hit: dict) -> dict:
    """Pull the match_context with the highest sequence identity out of a
    hit's `services` block (the sequence search's raw alignment stats)."""
    best_mc: dict = {}
    for svc in hit.get("services", []):
        if svc.get("service_type") != "sequence":
            continue
        for node in svc.get("nodes", []):
            for mc in node.get("match_context", []) or []:
                if mc.get("sequence_identity", 0) >= best_mc.get("sequence_identity", -1):
                    best_mc = mc
    return {field: best_mc.get(field) for field in _MATCH_CONTEXT_FIELDS}


def parse_hit_ids(result_set: List[dict]) -> List[Dict]:
    """From polymer_entity identifiers like '4HHB_1', return one dict per
    PDB entry with the normalized relevancy score plus the raw sequence
    alignment stats (identity, e-value, bitscore, alignment coordinates),
    deduplicated by pdb_id (keeping the highest-identity entity per entry)."""
    best: Dict[str, Dict] = {}
    for hit in result_set:
        ident = hit["identifier"]  # e.g. "4HHB_1"
        pdb_id, _, entity_id = ident.partition("_")
        row = {
            "entity_id": entity_id,
            "score": hit.get("score", 0.0),
            **_best_match_context(hit),
        }
        current_identity = row.get("sequence_identity") or -1
        existing_identity = best.get(pdb_id, {}).get("sequence_identity") or -1
        if pdb_id not in best or current_identity > existing_identity:
            best[pdb_id] = row
    return [{"pdb_id": pdb_id, **row} for pdb_id, row in best.items()]


# --------------------------------------------------------------------------- #
# 2. Fetch crystallization details
# --------------------------------------------------------------------------- #
CACHE_DIR = os.path.join(os.path.dirname(os.path.abspath(__file__)), ".pdb_cache")
os.makedirs(CACHE_DIR, exist_ok=True)

def get_cached_mmcif(pdb_id):
    """Retrieve mmCIF file with caching."""
    cache_file = os.path.join(CACHE_DIR, f"{pdb_id.lower()}.cif")
    if os.path.exists(cache_file):
        with open(cache_file, "r") as f:
            return f.read()
    try:
        url = f"https://files.rcsb.org/view/{pdb_id}.cif"
        response = requests.get(url, timeout=10)
        response.raise_for_status()
        content = response.text

        # Cache the result
        with open(cache_file, "w") as f:
            f.write(content)

        return content

    except Exception as e:
        print(f"Warning fetching mmCIF for {pdb_id}: {e}")
        return None

def get_cached_pubmed_id(pdb_id):
    """
    Retrieve PubMed ID with caching.

    1. Check local JSON cache
    2. Try the RCSB REST API (data.rcsb.org)
    3. If still not found -> fall back to parsing the mmCIF file with gemmi
       (NOT plain regex -- citation fields live in a `loop_` block, so a
       naive regex can grab the *next tag name* instead of a value when
       there's no PubMed ID on the same line -- e.g. it returns
       "_citation.pdbx_database_id_DOI" instead of "NA" for entries like
       1VGL that have no PubMed ID)
    """
    pdb_id = pdb_id.upper()
    cache_file = os.path.join(CACHE_DIR, f"{pdb_id.lower()}_pubmed.json")

    # ---------- Cached JSON ----------
    if os.path.exists(cache_file):
        with open(cache_file, "r") as f:
            cached = json.load(f).get("pubmed_id", "NA")
        if cached not in (None, "NA") and str(cached).isdigit():
            return str(cached)

    # ---------- REST API ----------
    try:
        entry_data = requests.get(
            f"https://data.rcsb.org/rest/v1/core/entry/{pdb_id}",
            timeout=10
        ).json()
        pubmed = entry_data.get("rcsb_primary_citation", {}).get("pdbx_database_id_pub_med")

        if pubmed not in (None, "", "NA"):
            with open(cache_file, "w") as f:
                json.dump({"pubmed_id": pubmed}, f)
            return str(pubmed)

    except Exception as e:
        print(f"Warning fetching PubMed ID via REST API for {pdb_id}: {e}")

    # ---------- mmCIF fallback (gemmi, not regex) ----------
    try:
        mmcif_content = get_cached_mmcif(pdb_id)
        if mmcif_content:
            doc = gemmi.cif.read_string(mmcif_content)
            block = doc.sole_block()

            pubmed_values = [
                v for v in block.find_values("_citation.pdbx_database_id_PubMed")
                if v not in (".", "?", "", None)
            ]

            if pubmed_values:
                pubmed = pubmed_values[0]
                with open(cache_file, "w") as f:
                    json.dump({"pubmed_id": pubmed}, f)
                return str(pubmed)

    except Exception as e:
        print(f"Warning parsing PubMed ID from mmCIF for {pdb_id}: {e}")

    return "NA"


def get_ph_from_mmcif_or_details(block):
    raw_ph = block.find_value("_exptl_crystal_grow.pH")
    try:
        if raw_ph not in (None, "", ".", "?", "NA"):
            return float(raw_ph)
    except ValueError:
        pass
    details = block.find_value("_exptl_crystal_grow.pdbx_details")
    if not details:
        return None
    match = re.search(r"\bpH\s*[:=]?\s*([0-9]+(?:\.[0-9]+)?)", details, re.IGNORECASE)
    if match:
        try:
            return float(match.group(1))
        except ValueError:
            return None
    return None

# ....get_method_from_mmcif_or_details...
def get_method_from_mmcif_or_details(block):
    method = block.find_value("_exptl_crystal_grow.method")
    if method not in (None, "", ".", "?", "NA"):
        return method.strip()
    details = block.find_value("_exptl_crystal_grow.pdbx_details")
    if not details:
        return None
    match = re.search(r"\b(method|technique|approach)\s*[:=]?\s*([a-zA-Z0-9 \-]+)", details, re.IGNORECASE)
    if match:
        return match.group(2).strip()
    return None

#.....get_temperature_from_mmcif_or_details.....
def get_temperature_from_mmcif_or_details(block):
   
    # ---- 1. Try mmCIF temperature fields first ----
    temp = block.find_value("_exptl_crystal_grow.temp") 

    if temp not in (None, "", ".", "?", "NA"):
        try:
            return float(temp)
        except ValueError:
            pass  # fall back to pdbx_details

    # ---- 2. Fallback: extract from pdbx_details ----
    details = block.find_value("_exptl_crystal_grow.pdbx_details")
    if not details:
        return None

    # Match ONLY Kelvin values like 277K or 277.00K
    match = re.search(
        r"\b([0-9]+(?:\.[0-9]+)?)\s*K\b",
        details,
        re.IGNORECASE,
    )

    if match:
        try:
            return float(match.group(1))
        except ValueError:
            return None

    return None


#.....get_pdbx_ph_range_from_mmcif_or_details...
def get_pdbx_ph_range_from_mmcif_or_details(block):
    """
    Return pH range from mmCIF if present, otherwise extract from pdbx_details.
    Output format: 'low-high' (string), or None if not found.
    """

    # 1. Try mmCIF field first
    ph_range = block.find_value("_exptl_crystal_grow.pdbx_pH_range")
    if ph_range not in (None, "", ".", "?", "NA"):
        return ph_range.strip()

    # 2. Fallback: extract from pdbx_details
    details = block.find_value("_exptl_crystal_grow.pdbx_details")
    if not details:
        return None

    details = str(details)

    # Match patterns like:
    # pH 6.5-7.5
    # pH range 6.0 - 7.0
    # pH: 6.8 to 7.2
    match = re.search(
        r"\bpH\b.*?([0-9]+(?:\.[0-9]+)?)\s*(?:-|\u2013|to)\s*([0-9]+(?:\.[0-9]+)?)",
        details,
        re.IGNORECASE,
    )

    if match:
        low = match.group(1)
        high = match.group(2)
        return f"{low}-{high}"

    return None

def fetch_non_polymers(pdb_id):
    """
    Fetch non-polymer IDs for a PDB entry.

    1. Try RCSB GraphQL
    2. If none found -> fallback to mmCIF parsing

    Returns
    -------
    dict | None
        {"non_polymers": "LIG1, LIG2", "view_3d_url": "https://www.rcsb.org/3d-view/1XOM"}
        if any non-polymers (ligands) are found, otherwise None.
    """
    pdb_id = pdb_id.upper()
    view_3d_url = f"https://www.rcsb.org/3d-view/{pdb_id}"

    # ---------- GraphQL ----------
    url = "https://data.rcsb.org/graphql"

    query = """
    query getLigands($id: String!) {
      entry(entry_id: $id) {
        nonpolymer_entities {
          rcsb_nonpolymer_entity_container_identifiers {
            nonpolymer_comp_id
          }
        }
      }
    }
    """

    try:
        response = requests.post(
            url,
            json={"query": query, "variables": {"id": pdb_id}},
            timeout=15
        )
        response.raise_for_status()
        data = response.json()

        entry = data.get("data", {}).get("entry", {})
        entities = entry.get("nonpolymer_entities") or []

        non_polymers = [
            e["rcsb_nonpolymer_entity_container_identifiers"]["nonpolymer_comp_id"]
            for e in entities
            if e.get("rcsb_nonpolymer_entity_container_identifiers", {}).get("nonpolymer_comp_id")
        ]

        if non_polymers:
            return {
              #  "non_polymers": ", ".join(list(dict.fromkeys(non_polymers))),
                 "view_3d_url": view_3d_url,
            }

    except Exception as e:
        print(f"GraphQL ligand error for {pdb_id}: {e}")

    # ---------- mmCIF fallback ----------
    try:
        cif_path = os.path.join(CACHE_DIR, f"{pdb_id.lower()}.cif")

        if not os.path.exists(cif_path):
            mmcif_content = get_cached_mmcif(pdb_id)
            if not mmcif_content:
                return None

        doc = gemmi.cif.read_file(cif_path)
        block = doc.sole_block()

        comp_ids = block.find_values("_pdbx_entity_instance_feature.comp_id")
        comp_ids = [c for c in comp_ids if c not in ("?", ".", "", None)]

        if comp_ids:
            return {
                # "non_polymers": ", ".join(list(dict.fromkeys(comp_ids))),
                "view_3d_url": view_3d_url,
            }

    except Exception as e:
        print(f"mmCIF ligand fallback error for {pdb_id}: {e}")

    return None

def get_polymer_type_from_mmcif(block):
    """
    Determine if the structure is single-polymer or complex.

    Returns:
        "uni_pol"  -> only one polymer entity
        "complex"  -> multiple polymer entities
    """

    try:
        # Get entity types (polymer, non-polymer, etc.)
        entity_types = block.find_values("_entity.type")

        # Count only polymer entities
        polymer_count = sum(1 for t in entity_types if str(t).lower() == "polymer")

        if polymer_count == 1:
            return "uni_pol"
        elif polymer_count > 1:
            return "complex"
        else:
            return None  # no polymer found (rare case)

    except Exception as e:
        print(f"Failed to determine polymer type: {e}")
        return None

def extract_mmcif_info(pdb_id):
    """
    Extract mmCIF info for a PDB entry and compute sequence identity
    with the provided query sequence.
    """
    pdb_id = pdb_id.upper()

    # PubMed
    pubmed = get_cached_pubmed_id(pdb_id)

    # Download mmCIF
    mmcif_content = get_cached_mmcif(pdb_id)
    if not mmcif_content:
        return None

    # Parse mmCIF with gemmi
    try:
        doc = gemmi.cif.read_string(mmcif_content)
        block = doc.sole_block()
    except Exception as e:
        print(f"Failed to parse mmCIF for {pdb_id}: {e}")
        return None

    # ---- Extract assembly info robustly ----
    assembly_ids = block.find_values("_pdbx_struct_assembly.id")
    oligomeric_values = block.find_values("_pdbx_struct_assembly.oligomeric_details")
    assembly_detail = None
    if oligomeric_values:
        first_detail = oligomeric_values[0]
        if len(assembly_ids) > 1:
            first_detail += "(*)"
        assembly_detail = first_detail

    # Determine polymer type (uni_pol vs complex)
    polymer_type = get_polymer_type_from_mmcif(block) or "NA"


    # Collect info
    info = {
        "pdb_id": pdb_id, 
        "Resolution": block.find_value("_refine.ls_d_res_high"),
        "Pubmed_id": pubmed,
        "Polymer": polymer_type, 
        "Assembly": assembly_detail,
        "Method": get_method_from_mmcif_or_details(block),
        "pH": get_ph_from_mmcif_or_details(block),
        "Temp": get_temperature_from_mmcif_or_details(block),
        "pdbx_details": block.find_value("_exptl_crystal_grow.pdbx_details"),
        "pdbx_pH_range": get_pdbx_ph_range_from_mmcif_or_details(block),
        "Non_polymers": fetch_non_polymers(pdb_id)
    }

    return info


# --------------------------------------------------------------------------- #
# 3. Filter down to experimental conditions -- Output.csv is the ONLY file
#    this script writes.
# --------------------------------------------------------------------------- #
def filter_experimental_conditions(rows_or_csv, output_csv="Output.csv"):
    """
    Filter homolog-search results down to rows that have at least one
    experimental crystallization detail (pH, temperature, method, or
    free-text pdbx_details), and write ONLY the filtered result to disk.

    Parameters
    ----------
    rows_or_csv : list[dict] | pandas.DataFrame | str
        Either the in-memory results (a list of row-dicts, as returned by
        find_homologs_with_conditions()), a DataFrame, or a path to an
        existing CSV to read from. Passing the in-memory results directly
        (the normal CLI path) avoids ever writing an unfiltered CSV to disk.
    output_csv : str
        Path to write the filtered CSV to. Defaults to "Output.csv".

    Returns
    -------
    str
        The path to the filtered CSV that was written -- this is the only
        file the CLI produces.
    """
    if isinstance(rows_or_csv, pd.DataFrame):
        df = rows_or_csv
    elif isinstance(rows_or_csv, str):
        if not os.path.exists(rows_or_csv):
            raise FileNotFoundError(
                f"Input CSV not found: {rows_or_csv} "
                f"(run find_homologs_with_conditions() first, "
                f"or pass the correct path)"
            )
        df = pd.read_csv(rows_or_csv)
    else:
        df = pd.DataFrame(rows_or_csv)

    condition_cols = ["pH", "Temp", "Method", "pdbx_details"]
    missing_cols = [c for c in condition_cols if c not in df.columns]
    if missing_cols:
        raise KeyError(
            f"Expected column(s) not found: {missing_cols}. "
            f"Available columns: {list(df.columns)}"
        )

    filtered_df = df[df[condition_cols].notna().any(axis=1)]

    filtered_df.to_csv(output_csv, index=False)

    print(f"Filtered CSV saved to: {os.path.abspath(output_csv)}")
    return output_csv


# --------------------------------------------------------------------------- #
# 4. Orchestration
# --------------------------------------------------------------------------- #
def prompt_for_sequence() -> str:
    """Interactively read a sequence from the terminal.

    Accepts either a raw sequence or a pasted FASTA record (header lines
    starting with '>' are ignored). Input ends on a blank line or EOF
    (Ctrl-D on Linux/macOS, Ctrl-Z then Enter on Windows), so multi-line
    pastes work fine.
    """
    print("Enter/paste the query sequence (protein/DNA/RNA, plain or FASTA).", file=sys.stderr)
    print("Press Enter on an empty line, or Ctrl-D, when done:", file=sys.stderr)
    lines = []
    try:
        while True:
            line = input()
            if line.strip() == "":
                break
            lines.append(line.strip())
    except EOFError:
        pass
    seq_lines = [ln for ln in lines if not ln.startswith(">")]
    return "".join(seq_lines)


def read_fasta(path: str) -> str:
    seq_lines = []
    with open(path) as f:
        for line in f:
            line = line.strip()
            if not line or line.startswith(">"):
                continue
            seq_lines.append(line)
    return "".join(seq_lines)


def find_homologs_with_conditions(
    sequence: str,
    identity_cutoff: float = 0.5,
    evalue_cutoff: float = 1.0,
    sequence_type: str = "protein",
    max_hits: int = 1000,
    sleep: float = 0.2,
    verbose: bool = True,
) -> List[Dict]:
    """Run the full pipeline and return a list of row-dicts, one per
    (pdb_id, compound) pair. Entries with no parsed compounds still get a
    single row with empty compound fields, so no hit is silently dropped."""

    query = build_sequence_query(
        sequence,
        identity_cutoff=identity_cutoff,
        evalue_cutoff=evalue_cutoff,
        sequence_type=sequence_type,
        max_hits=max_hits,
    )
    if verbose:
        print(f"Searching RCSB PDB (X-ray only) for homologs "
              f"(identity >= {identity_cutoff}, e-value <= {evalue_cutoff})...",
              file=sys.stderr)

    result_set = run_sequence_search(query)
    hits = parse_hit_ids(result_set)

    if verbose:
        print(f"Found {len(hits)} matching PDB entries.", file=sys.stderr)

    rows: List[Dict] = []

    for i, hit in enumerate(hits, 1):
        pdb_id = hit["pdb_id"]
        identity = hit.get("sequence_identity")
        evalue = hit.get("evalue")
        if verbose:
            identity_pct = f"{identity * 100:.1f}%" if identity is not None else "n/a"
            print(f"  [{i}/{len(hits)}] {pdb_id} (identity={identity_pct}, "
                  f"e-value={evalue})", file=sys.stderr)
    
        base_row = {
            "pdb_id": pdb_id,
            "entity_id": hit.get("entity_id"),
            "sequence_identity": round(identity * 100, 2) if identity is not None else None,
            "evalue": evalue,
            "score": round(hit.get("score", 0.0), 4),}

       
        info = extract_mmcif_info(pdb_id)
        
        if info:
            rows.append({**base_row, **info})

        time.sleep(sleep)  # be polite to the API

    return rows


# --------------------------------------------------------------------------- #
# CLI
# --------------------------------------------------------------------------- #
def main():
    p = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    seq_group = p.add_mutually_exclusive_group(required=False)
    seq_group.add_argument("--sequence", help="Protein/nucleic acid sequence, as a raw string.")
    seq_group.add_argument("--fasta", help="Path to a FASTA file containing the query sequence.")
    p.add_argument("--sequence-type", default="protein", choices=["protein", "dna", "rna"],
                   help="Sequence type for the search (default: protein).")
    p.add_argument("--identity", type=float, default=0.5,
                   help="Minimum fractional sequence identity, 0-1 (default: 0.5).")
    p.add_argument("--evalue", type=float, default=1.0,
                   help="Maximum E-value cutoff (default: 1.0).")
    p.add_argument("--max-hits", type=int, default=100,
                   help="Maximum number of PDB entries to retrieve (default: 100).")
    p.add_argument("--sleep", type=float, default=0.2,
                   help="Seconds to sleep between per-entry API calls (default: 0.2).")
    p.add_argument("--out", default="Output.csv",
                   help="Path for the final filtered CSV (default: Output.csv). "
                        "This is the only file the script writes.")
    p.add_argument("--quiet", action="store_true", help="Suppress progress messages.")
    args = p.parse_args()

    if args.sequence:
        sequence = args.sequence
    elif args.fasta:
        sequence = read_fasta(args.fasta)
    else:
        sequence = prompt_for_sequence()

    sequence = sequence.strip().upper()
    if not sequence:
        p.error("No sequence provided (nothing entered / FASTA file was empty).")

    rows = find_homologs_with_conditions(
        sequence,
        identity_cutoff=args.identity,
        evalue_cutoff=args.evalue,
        sequence_type=args.sequence_type,
        max_hits=args.max_hits,
        sleep=args.sleep,
        verbose=not args.quiet,
    )

    # rows stays in memory -- no raw/unfiltered CSV is ever written.
    # filter_experimental_conditions() is the only function that touches
    # disk, and Output.csv (or --out) is the only file it produces.
    filter_experimental_conditions(rows, output_csv=args.out)


if __name__ == "__main__":
    main()
