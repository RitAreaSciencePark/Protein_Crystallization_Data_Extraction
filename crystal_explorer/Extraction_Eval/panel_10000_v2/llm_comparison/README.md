# Open-weights LLM vs. the dictionary extractor

`llm_extract.py` runs an open-weights LLM on the 200 blinded audit entries
(`../second_annotator/`), or on all 1,000 with `--entries all`. The model
gets the same written guidelines as the human second annotator and returns a
JSON list of reagents. Decoding is greedy (temperature 0, seed 42).

`../second_annotator/score_agreement.py` then scores every finished run
exactly like the extractor: against annotator B and against the adjudicated
consensus. LLM items that cannot be matched automatically appear in
`adjudication.csv` for a human decision, like the extractor's.

Status: tested end to end against a mock server only. Ollama could not be
installed from the coding session (downloads were blocked there), and no
real model output exists yet.

## 1. Install Ollama (no admin rights needed)

```bash
mkdir -p ~/.local/ollama && cd ~/.local/ollama
curl -L -o ollama.tar.zst https://ollama.com/download/ollama-linux-amd64.tar.zst
tar --zstd -xf ollama.tar.zst && rm ollama.tar.zst
~/.local/ollama/bin/ollama serve &      # leave running
```

(With admin rights, `curl -fsSL https://ollama.com/install.sh | sh` installs
it as a service instead.)

## 2a. Pipeline test on the laptop (CPU, ~2 GB free RAM)

Only a ~1-2B model fits, and it is too weak to report. Use it only to check
that everything runs:

```bash
~/.local/ollama/bin/ollama pull qwen2.5:1.5b
cd crystal_explorer
.venv/bin/python3 Extraction_Eval/panel_10000_v2/llm_comparison/llm_extract.py --model qwen2.5:1.5b --limit 5
```

Delete `llm_runs/qwen2.5_1.5b/` afterwards if you do not want it scored.

## 2b. The real run, on a GPU node

Pick the strongest open model the GPU can hold, and check
ollama.com/library for current ones. As a rough guide, a ~70B model at 4-bit
needs about 48 GB of GPU memory, and an 8-14B model needs 16-24 GB.

```bash
ollama pull <model>          # e.g. llama3.3:70b or qwen2.5:72b
python llm_extract.py --model <model>                     # 200 entries
python llm_extract.py --model <model> --entries all       # optional: 1,000
```

With vLLM instead of Ollama (an OpenAI-compatible server):

```bash
vllm serve <hf-model-id> --port 8000
python llm_extract.py --api openai --base-url http://<node>:8000 --model <hf-model-id>
```

The script only needs `requests`. Copy this folder plus
`../second_annotator/{GUIDELINES.md,sample_meta.json}` and
`../audit_1000_raw.json` to the node, keeping the same relative layout, then
copy `llm_runs/<model>/` back here.

## 3. Score

After annotator B is done and A's list is reviewed (see
`../second_annotator/README.md`):

```bash
python ../second_annotator/score_agreement.py prepare   # adds the LLM's unmatched items
#   fill same_as / consensus for the new rows
python ../second_annotator/score_agreement.py score     # llm:<model>_vs_B, llm:<model>_vs_consensus
```

## Rules for a fair comparison

- **Don't edit the prompt after seeing scores on these entries.**
  `run_meta.json` stores its SHA-256 and `system_prompt.txt` stores the text.
  Tune on the development panel if needed.
- **Report everything about the run:** model name and size, quantization, the
  Ollama digest (in `run_meta.json`), the hardware, and the wall-clock time
  (`latency_s` per entry in `predictions.jsonl`).
- **Count failures as misses.** Entries that failed to parse are written with
  an `error` field and retried on the next run. A run is scored only when it
  covers all 200 entries; if some never succeed, record them as empty lists
  and say so in the paper.
