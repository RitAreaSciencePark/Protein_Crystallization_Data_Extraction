"""Run an open-weights LLM as a reagent extractor on the held-out audit
entries, for comparison with compound_extraction.py (paper Section VI).

The model gets exactly the written guidelines given to the human second
annotator (../second_annotator/GUIDELINES.md) and returns a JSON list of
reagents per entry. Decoding is greedy (temperature 0, fixed seed) so a run
can be repeated. Scoring is done by ../second_annotator/score_agreement.py,
which picks up every finished run in llm_runs/.

Backends
  --api ollama   Ollama's native /api/chat with a JSON-schema `format`
                 (local laptop test, or a GPU node running `ollama serve`)
  --api openai   any OpenAI-compatible server, e.g. vLLM
                 (`vllm serve <model>`), using `response_format` json_schema

Examples
  python llm_extract.py --model qwen2.5:1.5b --limit 5          # quick test
  python llm_extract.py --model llama3.3:70b                     # GPU node, 200 entries
  python llm_extract.py --api openai --base-url http://gpu01:8000 \\
         --model Qwen/Qwen2.5-72B-Instruct --entries all            # all 1,000

Held-out discipline: the prompt (guidelines + PROMPT_SUFFIX below) must not
be changed after looking at scores on these entries. Its SHA-256 is written
to run_meta.json; tune, if at all, on the development panel instead.

Resumable: entries already in predictions.jsonl are skipped.
"""
import argparse
import datetime
import hashlib
import json
import os
import re
import time

import requests

HERE = os.path.dirname(os.path.abspath(__file__))
PANEL = os.path.dirname(HERE)
GUIDELINES = os.path.join(PANEL, "second_annotator", "GUIDELINES.md")
SAMPLE_META = os.path.join(PANEL, "second_annotator", "sample_meta.json")
AUDIT = os.path.join(PANEL, "audit_1000_raw.json")

PROMPT_SUFFIX = """

---

You are the annotator. You will receive one `pdbx_details` text at a time.
Apply the guidelines above and answer with JSON only, in the form
{"reagents": ["name 1", "name 2", ...]}. Use an empty list if the entry
contains no reagent (the guidelines' `NONE`). Do not add explanations.
"""
SCHEMA = {"type": "object",
          "properties": {"reagents": {"type": "array", "items": {"type": "string"}}},
          "required": ["reagents"], "additionalProperties": False}
SEED = 42


def system_prompt():
    with open(GUIDELINES, encoding="utf-8") as f:
        text = f.read()
    # the human-facing "how to fill in the CSV" section is irrelevant to the model
    text = re.sub(r"## What to write.*?(?=## What counts)", "", text, flags=re.S)
    text = re.sub(r"When you are done, send back.*", "", text, flags=re.S)
    return text.strip() + PROMPT_SUFFIX


def load_entries(which):
    with open(AUDIT) as f:
        audit = {e["pdb_id"]: e for e in json.load(f)}
    if which == "all":
        ids = list(audit)
    else:
        with open(SAMPLE_META) as f:
            ids = json.load(f)["pdb_ids"]
    return [audit[i] for i in ids]


def call_ollama(base, model, system, text, timeout):
    r = requests.post(f"{base}/api/chat", timeout=timeout, json={
        "model": model, "stream": False, "format": SCHEMA,
        "options": {"temperature": 0, "seed": SEED, "num_ctx": 8192},
        "messages": [{"role": "system", "content": system}, {"role": "user", "content": text}],
    })
    r.raise_for_status()
    d = r.json()
    return d["message"]["content"], d.get("prompt_eval_count"), d.get("eval_count")


def call_openai(base, model, system, text, timeout):
    r = requests.post(f"{base}/v1/chat/completions", timeout=timeout, json={
        "model": model, "temperature": 0, "seed": SEED,
        "response_format": {"type": "json_schema",
                            "json_schema": {"name": "reagents", "schema": SCHEMA, "strict": True}},
        "messages": [{"role": "system", "content": system}, {"role": "user", "content": text}],
    })
    r.raise_for_status()
    d = r.json()
    u = d.get("usage") or {}
    return d["choices"][0]["message"]["content"], u.get("prompt_tokens"), u.get("completion_tokens")


def parse(raw):
    try:
        items = json.loads(raw)["reagents"]
    except (json.JSONDecodeError, KeyError, TypeError):
        m = re.search(r"\{.*\}", raw or "", re.S)  # tolerate stray text around the JSON
        items = json.loads(m.group(0))["reagents"] if m else None
    if not isinstance(items, list):
        raise ValueError("no reagent list in model output")
    out = []
    for x in items:
        x = str(x).strip()
        if x and x.upper() != "NONE" and x not in out:
            out.append(x)
    return out


def model_info(api, base, model):
    if api != "ollama":
        return {}
    try:
        d = requests.post(f"{base}/api/show", json={"model": model}, timeout=30).json()
        det = d.get("details", {})
        return {"family": det.get("family"), "parameter_size": det.get("parameter_size"),
                "quantization": det.get("quantization_level"),
                "digest": next((m.get("digest") for m in requests.get(f"{base}/api/tags", timeout=30)
                                .json().get("models", []) if m.get("name") == model), None)}
    except requests.RequestException:
        return {}


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--model", required=True)
    ap.add_argument("--api", choices=["ollama", "openai"], default="ollama")
    ap.add_argument("--base-url", default="http://localhost:11434")
    ap.add_argument("--entries", choices=["sample200", "all"], default="sample200")
    ap.add_argument("--limit", type=int, default=None, help="only the first N entries (testing)")
    ap.add_argument("--timeout", type=float, default=600)
    args = ap.parse_args()

    system = system_prompt()
    entries = load_entries(args.entries)[: args.limit]
    run_dir = os.path.join(HERE, "llm_runs", re.sub(r"[^A-Za-z0-9._-]", "_", args.model))
    os.makedirs(run_dir, exist_ok=True)
    pred_path = os.path.join(run_dir, "predictions.jsonl")
    done = set()
    if os.path.exists(pred_path):
        with open(pred_path) as f:
            done = {json.loads(l)["pdb_id"] for l in f if l.strip() and '"error"' not in l}

    meta = {"model": args.model, "api": args.api, "base_url": args.base_url, "entries": args.entries,
            "temperature": 0, "seed": SEED, "prompt_sha256": hashlib.sha256(system.encode()).hexdigest(),
            "started": datetime.datetime.now().isoformat(timespec="seconds"),
            **model_info(args.api, args.base_url, args.model)}
    with open(os.path.join(run_dir, "run_meta.json"), "w") as f:
        json.dump(meta, f, indent=2)
    with open(os.path.join(run_dir, "system_prompt.txt"), "w", encoding="utf-8") as f:
        f.write(system)

    call = call_ollama if args.api == "ollama" else call_openai
    todo = [e for e in entries if e["pdb_id"] not in done]
    print(f"{args.model}: {len(entries)} entries, {len(todo)} to do", flush=True)
    t_all = time.time()
    with open(pred_path, "a") as out:
        for k, e in enumerate(todo, 1):
            text = e["pdbx_details"].strip("'").strip(";").strip()
            t0 = time.time()
            rec = {"pdb_id": e["pdb_id"]}
            try:
                raw, p_tok, c_tok = call(args.base_url, args.model, system, text, args.timeout)
                rec.update(reagents=parse(raw), raw=raw, prompt_tokens=p_tok, completion_tokens=c_tok)
            except Exception as ex:  # one bad entry must not stop the run; it is retried next time
                rec.update(error=f"{type(ex).__name__}: {ex}")
            rec["latency_s"] = round(time.time() - t0, 2)
            out.write(json.dumps(rec) + "\n")
            out.flush()
            print(f"[{k}/{len(todo)}] {e['pdb_id']} {rec['latency_s']}s "
                  f"{rec.get('reagents', rec.get('error'))}", flush=True)
    print(f"done in {time.time() - t_all:.0f}s -> {pred_path}")


if __name__ == "__main__":
    main()
