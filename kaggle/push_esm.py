#!/usr/bin/env python3
"""Build and push the private ESM scoring kernel, or fetch its output.

    python kaggle/push_esm.py push              # embed the scoring input, push to Kaggle
    python kaggle/push_esm.py status
    python kaggle/push_esm.py fetch --dest data/plm/esm_scores

Uses the caller's own kaggle CLI configuration. The kernel is private, runs on
one T4 GPU with internet (to download the pinned Hugging Face revisions) and
receives only protein sequences and positions: no labels, no variant IDs.
"""
from __future__ import annotations

import argparse
import base64
import gzip
import hashlib
import json
import os
import subprocess
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
INPUT = ROOT / "data/plm/esm_scoring_input.json"
KERNEL = ROOT / "kaggle/esm_score_kernel.py"
SLUG = "hl-esm-masked-marginal"
BUILD = Path(os.environ.get("KAGGLE_BUILD_ROOT", Path.home() / "work/hlplm/kaggle_build")) / SLUG


def kaggle(*args: str) -> str:
    print("$ kaggle", *args, flush=True)
    return subprocess.run(["kaggle", *args], check=True, text=True, capture_output=True).stdout


def username() -> str:
    for line in kaggle("config", "view").splitlines():
        if line.strip().startswith("- username:"):
            return line.split(":", 1)[1].strip()
    raise SystemExit("kaggle CLI has no configured username")


def build(user: str) -> Path:
    raw = INPUT.read_bytes()
    text = KERNEL.read_text(encoding="utf-8")
    if text.count("__INPUT_B64__") != 1 or text.count("__INPUT_SHA256__") != 1:
        raise SystemExit("kernel placeholders missing")
    payload = base64.b64encode(gzip.compress(raw, mtime=0)).decode()
    text = (text.replace("__INPUT_B64__", payload)
            .replace("__INPUT_SHA256__", hashlib.sha256(raw).hexdigest()))
    BUILD.mkdir(parents=True, exist_ok=True)
    (BUILD / KERNEL.name).write_text(text, encoding="utf-8")
    meta = {"id": f"{user}/{SLUG}", "title": SLUG, "code_file": KERNEL.name, "language": "python",
            "kernel_type": "script", "is_private": True, "enable_gpu": True,
            "enable_internet": True, "machine_shape": "NvidiaTeslaT4", "dataset_sources": [],
            "kernel_sources": [], "competition_sources": []}
    (BUILD / "kernel-metadata.json").write_text(json.dumps(meta, indent=2) + "\n")
    print(f"built {BUILD} (input sha256 {hashlib.sha256(raw).hexdigest()})")
    return BUILD


def main() -> None:
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("action", choices=["build", "push", "status", "fetch"])
    ap.add_argument("--dest", type=Path, default=ROOT / "data/plm/esm_scores")
    args = ap.parse_args()
    user = username()
    if args.action in ("build", "push"):
        path = build(user)
        if args.action == "push":
            print(kaggle("kernels", "push", "-p", str(path)))
    elif args.action == "status":
        print(kaggle("kernels", "status", f"{user}/{SLUG}"))
    else:
        args.dest.mkdir(parents=True, exist_ok=True)
        print(kaggle("kernels", "output", f"{user}/{SLUG}", "-p", str(args.dest), "-o"))


if __name__ == "__main__":
    main()
