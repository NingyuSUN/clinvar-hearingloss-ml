"""Masked-marginal ESM scoring for the PLM-ESM-2026-10 protocol (Kaggle GPU kernel).

``kaggle/push_esm.py`` embeds ``data/plm/esm_scoring_input.json`` below and
pushes this file as a private Kaggle script. For each requested model, every
(protein, position) pair is scored by masking that one residue inside its
prespecified window and recording log-softmax probabilities of the 20 standard
residues. No labels or variant identities are present in the input.

Outputs (in the working directory): ``esm_logprobs.csv.gz`` and
``esm_scoring_manifest.json`` (versions, revisions, weight/input/output SHA256).

Local smoke test with a tiny model (CPU):
    python kaggle/esm_score_kernel.py --input data/plm/esm_scoring_input.json \
        --models facebook/esm2_t6_8M_UR50D --limit 20 --out ~/tmp/esm_smoke
"""
from __future__ import annotations

import argparse
import base64
import gc
import gzip
import hashlib
import io
import json
import os
import platform
import queue
import shutil
import threading
import time
from datetime import datetime, timezone
from pathlib import Path

AA20 = "ACDEFGHIKLMNPQRSTVWY"
INPUT_B64 = "__INPUT_B64__"
INPUT_SHA256 = "__INPUT_SHA256__"
MODEL_REVISIONS = {
    "facebook/esm1v_t33_650M_UR90S_1": "8bfdb1892536cc77bd0760b9c25ddced2cd0b4c8",
    "facebook/esm1v_t33_650M_UR90S_2": "3c1e9e64480f069b163e456c361bfd80a8bab04c",
    "facebook/esm1v_t33_650M_UR90S_3": "0b00fd112e63f6b5e70a9cd8484d4e660312ce70",
    "facebook/esm1v_t33_650M_UR90S_4": "443968f644da132d323bbc6321a6d443149fa57b",
    "facebook/esm1v_t33_650M_UR90S_5": "fb2e51cb0f605cbad2c4ca3bb784be7bc8e2f4a8",
    "facebook/esm2_t33_650M_UR50D": "08e4846e537177426273712802403f7ba8261b6c",
}
WEIGHT_FILES = ("model.safetensors", "pytorch_model.bin")
TOKENS_PER_BATCH = 8192
COLUMNS = ["transcript_stable_id", "protein_position", "model", "model_revision",
           "window_start", "window_end", "ref_aa", *[f"lp_{a}" for a in AA20]]


def sha256_bytes(data: bytes) -> str:
    return hashlib.sha256(data).hexdigest()


def sha256_file(path: Path) -> str:
    h = hashlib.sha256()
    with open(path, "rb") as fh:
        for chunk in iter(lambda: fh.read(1 << 22), b""):
            h.update(chunk)
    return h.hexdigest()


def load_input(path: str | None) -> tuple[dict, str]:
    if path:
        raw = Path(path).read_bytes()
    else:
        raw = gzip.decompress(base64.b64decode(INPUT_B64))
        if sha256_bytes(raw) != INPUT_SHA256:
            raise SystemExit("embedded scoring input does not match its SHA256")
    return json.loads(raw), sha256_bytes(raw)


def work_items(spec: dict, limit: int | None) -> list[dict]:
    items = []
    for prot in spec["proteins"]:
        seq = prot["sequence"]
        for pos in prot["positions"]:
            start, end = pos["window"]
            p = pos["position"]
            if not (0 <= start < p <= end <= len(seq)) or end - start > spec["max_residues"]:
                raise SystemExit(f"bad window for {prot['transcript_stable_id']}:{p}")
            items.append({"tid": prot["transcript_stable_id"], "position": p,
                          "start": start, "end": end, "window_seq": seq[start:end],
                          "ref": seq[p - 1]})
    items.sort(key=lambda r: (r["tid"], r["start"], r["position"]))
    return items[:limit] if limit else items


def fetch_model(repo: str, revision: str | None, cache: Path) -> tuple[Path, dict]:
    from huggingface_hub import HfApi, snapshot_download

    info = HfApi().model_info(repo, revision=revision, files_metadata=True)
    names = {s.rfilename for s in info.siblings}
    weight = next(w for w in WEIGHT_FILES if w in names)
    patterns = ["config.json", "vocab.txt", "special_tokens_map.json", "tokenizer_config.json", weight]
    local = Path(snapshot_download(repo, revision=info.sha, allow_patterns=patterns,
                                   cache_dir=str(cache)))
    meta = {"repo": repo, "requested_revision": revision, "resolved_revision": info.sha,
            "weight_file": weight, "weight_sha256": sha256_file(local / weight),
            "files": sorted(p.name for p in local.iterdir())}
    lfs = next((s.lfs for s in info.siblings if s.rfilename == weight), None)
    if lfs is not None and getattr(lfs, "sha256", None):
        meta["hub_lfs_sha256"] = lfs.sha256
        if lfs.sha256 != meta["weight_sha256"]:
            raise SystemExit(f"{repo}: downloaded weights do not match the hub LFS SHA256")
    return local, meta


def score_model(model_dir: Path, items: list[dict], device: str) -> list[dict]:
    import torch
    from transformers import EsmForMaskedLM, EsmTokenizer

    tok = EsmTokenizer.from_pretrained(str(model_dir))
    model = EsmForMaskedLM.from_pretrained(str(model_dir), torch_dtype=torch.float32).to(device).eval()
    aa_ids = tok.convert_tokens_to_ids(list(AA20))
    unk = tok.unk_token_id
    if unk in aa_ids:
        raise SystemExit("tokenizer lacks a standard residue")
    rows = []
    # Group positions sharing a window so one token template serves a batch.
    groups: dict[tuple, list[dict]] = {}
    for it in items:
        groups.setdefault((it["tid"], it["start"], it["end"]), []).append(it)
    with torch.inference_mode():
        for (tid, start, end), members in groups.items():
            window = members[0]["window_seq"]
            ids = [tok.cls_token_id, *tok.convert_tokens_to_ids(list(window)), tok.eos_token_id]
            base = torch.tensor(ids, dtype=torch.long)
            per_batch = max(1, TOKENS_PER_BATCH // len(ids))
            for i in range(0, len(members), per_batch):
                chunk = members[i:i + per_batch]
                batch = base.repeat(len(chunk), 1)
                idx = torch.tensor([m["position"] - start for m in chunk])  # +1 for <cls>, -1 for 0-based
                if not all(window[m["position"] - start - 1] == m["ref"] for m in chunk):
                    raise SystemExit(f"window residue mismatch in {tid}")
                batch[torch.arange(len(chunk)), idx] = tok.mask_token_id
                logits = model(input_ids=batch.to(device),
                               attention_mask=torch.ones_like(batch).to(device)).logits
                lp = torch.log_softmax(logits[torch.arange(len(chunk)), idx.to(device)].float(), dim=-1)
                lp = lp[:, aa_ids].cpu().numpy()
                for m, row in zip(chunk, lp):
                    rows.append({"transcript_stable_id": tid, "protein_position": m["position"],
                                 "window_start": start, "window_end": end, "ref_aa": m["ref"],
                                 **{f"lp_{a}": float(v) for a, v in zip(AA20, row)}})
    return rows


def release(device: str) -> None:
    """Free a finished model. HF modules hold reference cycles, so ``del`` alone
    left five 650M models resident and the sixth load ran out of GPU memory."""
    import torch

    gc.collect()
    if device.startswith("cuda"):
        with torch.cuda.device(device):
            torch.cuda.empty_cache()


def write_rows(path: Path, rows: list[dict]) -> None:
    lines = [",".join(COLUMNS)]
    for r in sorted(rows, key=lambda r: (r["model"], r["transcript_stable_id"], r["protein_position"])):
        lines.append(",".join(f"{r[c]:.6f}" if isinstance(r[c], float) else str(r[c]) for c in COLUMNS))
    buf = io.BytesIO()
    with gzip.GzipFile(fileobj=buf, mode="wb", mtime=0) as gz:
        gz.write(("\n".join(lines) + "\n").encode())
    path.write_bytes(buf.getvalue())


def run_model(repo: str, items: list[dict], device: str, cache_root: Path, out: Path):
    t0 = time.time()
    cache = cache_root / repo.replace("/", "__")
    local, meta = fetch_model(repo, MODEL_REVISIONS.get(repo), cache)
    rows = score_model(local, items, device)
    release(device)
    shutil.rmtree(cache, ignore_errors=True)  # keep disk use to the models in flight
    for r in rows:
        r["model"] = repo
        r["model_revision"] = meta["resolved_revision"]
    write_rows(out / f"partial__{repo.split('/')[-1]}.csv.gz", rows)  # survives a later failure
    meta |= {"n_positions": len(rows), "seconds": round(time.time() - t0, 1), "device": device}
    print(json.dumps(meta), flush=True)
    return rows, meta


def main() -> None:
    ap = argparse.ArgumentParser()
    ap.add_argument("--input", default=None, help="scoring input JSON (default: embedded copy)")
    ap.add_argument("--models", nargs="*", default=list(MODEL_REVISIONS))
    ap.add_argument("--limit", type=int, default=None, help="smoke test: first N positions")
    ap.add_argument("--out", default=".")
    ap.add_argument("--cache", default=os.environ.get("ESM_CACHE", "/tmp/hf_esm"))
    args = ap.parse_args()

    import numpy as np
    import torch
    import transformers
    import huggingface_hub

    torch.manual_seed(0)
    torch.backends.cuda.matmul.allow_tf32 = False
    torch.backends.cudnn.allow_tf32 = False
    n_gpu = torch.cuda.device_count()
    devices = [f"cuda:{i}" for i in range(n_gpu)] or ["cpu"]
    out = Path(args.out).expanduser()
    out.mkdir(parents=True, exist_ok=True)
    spec, input_sha = load_input(args.input)
    items = work_items(spec, args.limit)
    print(f"{len(items)} positions on {len({i['tid'] for i in items})} proteins; devices={devices}",
          flush=True)

    # One worker per GPU pulls models from a queue; results do not depend on
    # which device scored a model, and the output is sorted afterwards.
    todo: queue.Queue = queue.Queue()
    for repo in args.models:
        todo.put(repo)
    results, errors = {}, []

    def worker(device: str) -> None:
        while not errors:
            try:
                repo = todo.get_nowait()
            except queue.Empty:
                return
            try:
                results[repo] = run_model(repo, items, device, Path(args.cache), out)
            except BaseException as exc:  # surfaced after join
                errors.append((repo, device, repr(exc)))

    threads = [threading.Thread(target=worker, args=(d,)) for d in devices]
    for t in threads:
        t.start()
    for t in threads:
        t.join()
    if errors or set(results) != set(args.models):
        raise SystemExit(f"scoring failed: {errors or sorted(set(args.models) - set(results))}")
    all_rows = [r for repo in args.models for r in results[repo][0]]
    models_meta = [results[repo][1] for repo in args.models]
    out_path = out / "esm_logprobs.csv.gz"
    write_rows(out_path, all_rows)

    manifest = {
        "protocol": spec.get("protocol"),
        "generated_utc": datetime.now(timezone.utc).isoformat(timespec="seconds"),
        "strategy": "masked_marginal: one masked residue per forward pass, fp32, log_softmax over full vocabulary",
        "input_sha256": input_sha, "embedded_input": args.input is None,
        "ensembl_release": spec.get("ensembl_release"), "max_residues": spec.get("max_residues"),
        "n_positions_requested": len(items), "limit": args.limit,
        "models": models_meta,
        "output_file": out_path.name, "output_sha256": sha256_file(out_path),
        "output_uncompressed_sha256": sha256_bytes(gzip.decompress(out_path.read_bytes())),
        "versions": {"python": platform.python_version(), "torch": torch.__version__,
                     "transformers": transformers.__version__,
                     "huggingface_hub": huggingface_hub.__version__, "numpy": np.__version__,
                     "cuda": torch.version.cuda,
                     "devices": [torch.cuda.get_device_name(i) for i in range(n_gpu)] or ["cpu"]},
    }
    (out / "esm_scoring_manifest.json").write_text(json.dumps(manifest, indent=2) + "\n")
    print(json.dumps({k: v for k, v in manifest.items() if k != "models"}, indent=2))


if __name__ == "__main__":
    main()
