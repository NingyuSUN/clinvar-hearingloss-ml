#!/usr/bin/env python3
"""Validate the committed public artifacts in this repository.

This is a fast, dependency-free (stdlib only) sanity check meant to run in CI
on every push, without retraining anything. It exists because of a concrete
incident: a Windows-side tool once re-saved every file under `results/` with
CRLF line endings, which git reported as 17 "modified" files even though the
numeric content was byte-identical. This script catches that class of problem
(and a few other cheap-to-check regressions) before it reaches `main`.

Usage:
    python tools/validate_public_artifacts.py
"""

from __future__ import annotations

import csv
import json
import re
import sys
from pathlib import Path

REPO_ROOT = Path(__file__).resolve().parent.parent

# Files/dirs whose line endings we care about (text artifacts that get
# re-generated or hand-edited and could silently pick up CRLF).
TEXT_GLOBS = [
    "README.md",
    "docs/**/*.md",
    "results/**/*.csv",
    "results/**/*.json",
    "data/README.md",
    "*.cff",
]

REQUIRED_SUMMARY_KEYS = {
    "pooled_headline",
    "pooled_missense",
    "split_gap",
    "review_status_sensitivity",
}

REQUIRED_MANIFEST_KEYS = {
    "started_utc",
    "protocol",
    "k_outer",
    "inner_valid_frac",
    "num_boost_round_cap",
    "model_selection",
    "versions",
    "dataset_manifest_sha256",
    "seconds",
    "tags",
}

SHA256_RE = re.compile(r"^[0-9a-f]{64}$")


def fail(problems: list[str], message: str) -> None:
    problems.append(message)


def check_line_endings(problems: list[str]) -> None:
    seen: set[Path] = set()
    for pattern in TEXT_GLOBS:
        for path in REPO_ROOT.glob(pattern):
            if not path.is_file() or path in seen:
                continue
            seen.add(path)
            raw = path.read_bytes()
            if b"\r\n" in raw:
                fail(
                    problems,
                    f"CRLF line endings in {path.relative_to(REPO_ROOT)} "
                    "(expected LF; check .gitattributes normalization)",
                )


def check_json_file(problems: list[str], relpath: str, required_keys: set[str]) -> dict | None:
    path = REPO_ROOT / relpath
    if not path.exists():
        fail(problems, f"missing expected artifact: {relpath}")
        return None
    try:
        data = json.loads(path.read_text(encoding="utf-8"))
    except json.JSONDecodeError as exc:
        fail(problems, f"{relpath} is not valid JSON: {exc}")
        return None
    missing = required_keys - data.keys()
    if missing:
        fail(problems, f"{relpath} is missing expected keys: {sorted(missing)}")
    return data


def check_manifest_hash(problems: list[str], manifest: dict | None) -> None:
    if not manifest:
        return
    digest = manifest.get("dataset_manifest_sha256", "")
    if not SHA256_RE.match(digest):
        fail(
            problems,
            "results/run_manifest.json: dataset_manifest_sha256 is not a "
            "64-character hex sha256 digest",
        )


def check_result_csvs(problems: list[str]) -> None:
    csv_paths = sorted((REPO_ROOT / "results").rglob("*.csv"))
    if not csv_paths:
        fail(problems, "no CSV files found under results/")
        return
    for path in csv_paths:
        with path.open(newline="", encoding="utf-8") as fh:
            reader = csv.reader(fh)
            try:
                header = next(reader)
            except StopIteration:
                fail(problems, f"{path.relative_to(REPO_ROOT)} is empty (no header row)")
                continue
            if not header or any(col.strip() == "" for col in header):
                fail(problems, f"{path.relative_to(REPO_ROOT)} has a blank column header")
            n_data_rows = sum(1 for row in reader if row)
            if n_data_rows == 0:
                fail(problems, f"{path.relative_to(REPO_ROOT)} has a header but no data rows")


def check_no_stray_bytecode(problems: list[str]) -> None:
    for path in REPO_ROOT.rglob("__pycache__"):
        if ".git" in path.parts:
            continue
        tracked_marker = REPO_ROOT / ".git"
        if tracked_marker.exists():
            import subprocess

            out = subprocess.run(
                ["git", "ls-files", "--error-unmatch", str(path)],
                cwd=REPO_ROOT,
                capture_output=True,
                text=True,
            )
            if out.returncode == 0:
                fail(problems, f"__pycache__ directory is tracked by git: {path}")


def main() -> int:
    problems: list[str] = []

    check_line_endings(problems)
    summary = check_json_file(problems, "results/summary.json", REQUIRED_SUMMARY_KEYS)
    manifest = check_json_file(problems, "results/run_manifest.json", REQUIRED_MANIFEST_KEYS)
    check_manifest_hash(problems, manifest)
    check_result_csvs(problems)
    check_no_stray_bytecode(problems)

    if problems:
        print(f"validate_public_artifacts: {len(problems)} problem(s) found:\n")
        for problem in problems:
            print(f"  - {problem}")
        return 1

    n_csv = len(list((REPO_ROOT / "results").rglob("*.csv")))
    print(
        "validate_public_artifacts: OK "
        f"(summary keys={sorted(summary) if summary else []}, "
        f"{n_csv} result CSVs, LF line endings, manifest hash well-formed)"
    )
    return 0


if __name__ == "__main__":
    sys.exit(main())
