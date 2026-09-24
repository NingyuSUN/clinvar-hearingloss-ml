"""Malformed/adversarial input fuzzing for the predict/ CLIs and parser.

Contract under test: a bad *row* must come back as that row's error status
and never abort the batch; a bad *file* (unreadable, missing or ambiguous
header) must be rejected up front with a non-zero exit and no partial output.

All fuzzing is seeded, so a failure reproduces exactly. Needs xgboost/numpy
(predict/requirements-inference.txt), same as test_predict.py.
"""

from __future__ import annotations

import csv
import io
import json
import random
import re
import sys
from pathlib import Path

import pytest

REPO_ROOT = Path(__file__).resolve().parent.parent
sys.path.insert(0, str(REPO_ROOT / "predict"))

pytest.importorskip("xgboost")
pytest.importorskip("numpy")

import prediction_core as pc  # noqa: E402
import predict_variants  # noqa: E402

BUNDLE_DIR = REPO_ROOT / "predict" / "bundle"
HEADER = ",".join(pc.INPUT_FIELDS)
KNOWN_ROW = "48,GRCh38,16,2496587,G,C"  # a scorable frozen-cache variant

NORMALIZE_ERRORS = {
    "missing_variant_id",
    "unsupported_assembly",
    "unsupported_chromosome",
    "invalid_position",
    "invalid_allele",
    "identical_alleles",
    "unsupported_non_snv",
}
KEY_PATTERN = re.compile(r"(?:[1-9]|1[0-9]|2[0-2]|X|Y|MT)-[1-9][0-9]{0,9}-[ACGT]+-[ACGT]+")

# Hand-picked hostile values; the generator below mixes these with random junk.
NASTY = [
    None, "", " ", "\t", "\n", "\r\n", "\x00", "﻿", "​", " ",
    "0", "-1", "+5", "1.5", "1e9", "0x10", "00012", "NaN", "inf", "１２３", "٣",
    "1" * 11, "9" * 5000, "chr", "chrM", "chrUn", "M", "MT", "x", "23", "0",
    "A", "a", "N", "U", "ACGTN", "-", ".", "<DEL>", "A" * 1000, "é", "Ω", "💥",
    "GRCh38 ", " grch38", "GRCH38", "hg38", "GRCh37", "'; DROP TABLE variants; --",
    "${jndi:x}", "{{7*7}}", "../../etc/passwd", "\\", '"', "','", "\x1b[31m",
    5, 0, -3, 1.0, True, [], {}, ["1"],
]


def _junk(rng):
    kind = rng.random()
    if kind < 0.5:
        return rng.choice(NASTY)
    if kind < 0.75:
        return "".join(chr(rng.randrange(0, 0x3000)) for _ in range(rng.randrange(0, 12)))
    return "".join(rng.choice("ACGTacgtN0123456789chrXYM-") for _ in range(rng.randrange(0, 12)))


def _fuzz_row(rng, index):
    """Start from a valid-looking row and corrupt a random subset of fields."""
    row = {"variant_id": f"fz{index}", "assembly": "GRCh38", "chrom": "1",
           "pos": str(900_000_000 + index), "ref": "A", "alt": "G"}
    for field in rng.sample(pc.INPUT_FIELDS, rng.randrange(1, len(pc.INPUT_FIELDS) + 1)):
        if rng.random() < 0.1:
            del row[field]
        else:
            row[field] = _junk(rng)
    return row


def test_normalize_variant_never_raises_and_only_emits_canonical_keys():
    rng = random.Random(20260924)
    for i in range(20_000):
        row = _fuzz_row(rng, i)
        result = pc.normalize_variant(row)
        assert result["error"] is None or result["error"] in NORMALIZE_ERRORS, (row, result)
        if result["variant_key"] is not None:
            assert KEY_PATTERN.fullmatch(result["variant_key"]), (row, result)
            assert isinstance(result["pos"], int) and result["pos"] > 0
        if result["error"] is None:
            assert result["variant_key"] is not None and len(result["ref"]) == len(result["alt"]) == 1


def test_oversized_position_is_a_row_error_not_a_crash():
    row = {"variant_id": "p", "assembly": "GRCh38", "chrom": "1", "pos": "1" * 5000, "ref": "A", "alt": "G"}
    assert pc.normalize_variant(row)["error"] == "invalid_position"


@pytest.fixture(scope="module")
def bundle():
    return pc.Bundle(BUNDLE_DIR)


def test_fuzzed_batch_keeps_every_row_and_does_not_contaminate_good_rows(bundle):
    rng = random.Random(7)
    fuzzed = [_fuzz_row(rng, i) for i in range(5_000)]
    good_id, *good_fields = KNOWN_ROW.split(",")
    good = dict(zip(pc.INPUT_FIELDS, ["known-row", *good_fields]))
    rows = fuzzed[:2_500] + [good] + fuzzed[2_500:]

    solo = pc.predict_batch([good], bundle)[0]
    results = pc.predict_batch(rows, bundle)

    assert len(results) == len(rows)
    assert [r["input_index"] for r in results] == list(range(len(rows)))
    for r in results:
        assert r["annotation_status"], r
        for mo in r["models"].values():
            assert mo["status"] != "not_scored", r
            assert mo["probability"] is None
    buried = results[2_500]
    assert buried["variant_id"] == "known-row"
    for name in pc.MODELS:
        assert buried["models"][name]["status"] == "scored_research"
        assert buried["models"][name]["score"] == solo["models"][name]["score"]


def _run(cli, input_bytes, tmp_path, monkeypatch, bundle):
    """Invoke a CLI's main() in-process, reusing the already-verified bundle."""
    input_path = tmp_path / "input.csv"
    input_path.write_bytes(input_bytes)
    out = tmp_path / "out"
    monkeypatch.setattr(sys, "argv", [cli.__name__, "--input", str(input_path),
                                      "--bundle", str(BUNDLE_DIR), "--output", str(out)])
    monkeypatch.setattr(cli, "Bundle", lambda path: bundle)
    monkeypatch.delenv("ALPHAGENOME_API_KEY", raising=False)
    cli.main()
    return out


def _expected_rows(input_bytes):
    return list(csv.DictReader(io.StringIO(input_bytes.decode("utf-8-sig"), newline="")))


ROW_TOLERANT_FILES = {
    "bom": ("﻿" + HEADER + "\n" + KNOWN_ROW + "\n").encode(),
    "crlf": (HEADER + "\r\n" + KNOWN_ROW + "\r\n").encode(),
    "blank_lines": (HEADER + "\n\n" + KNOWN_ROW + "\n\n\n").encode(),
    "ragged_short_rows": (HEADER + "\n" + KNOWN_ROW + "\nshort,GRCh38\nonly_id\n").encode(),
    "extra_fields": (HEADER + "\n" + KNOWN_ROW + ",extra,more\n").encode(),
    "extra_named_column": (HEADER + ",notes\n" + KNOWN_ROW + ",hello\n").encode(),
    "reordered_columns": ("alt,ref,pos,chrom,assembly,variant_id\nC,G,2496587,16,GRCh38,48\n").encode(),
    "quoted_newline_in_id": (HEADER + "\n" + KNOWN_ROW + '\n"a\nb",GRCh38,1,5,A,G\n').encode(),
    "nul_bytes": (HEADER + "\n" + KNOWN_ROW + "\nx\x00,GRCh38,1\x00,5,A,G\n").encode(),
    "unterminated_quote": (HEADER + "\n" + KNOWN_ROW + '\n"open,GRCh38,1,5,A,G\n').encode(),
    "oversized_position": (HEADER + "\n" + KNOWN_ROW + "\np,GRCh38,1," + "1" * 5000 + ",A,G\n").encode(),
    "whitespace_padding": (HEADER + "\n 48 , GRCh38 , chr16 , 2496587 , g , c \n").encode(),
    "injection_strings": (HEADER + "\n" + KNOWN_ROW + '\n"=cmd|\'/c calc\'!A1",GRCh38,1,5,A,G\n').encode(),
}


@pytest.mark.parametrize("name", sorted(ROW_TOLERANT_FILES))
def test_cli_row_level_malformations_are_per_row(name, tmp_path, monkeypatch, bundle):
    data = ROW_TOLERANT_FILES[name]
    out = _run(predict_variants, data, tmp_path, monkeypatch, bundle)

    records = [json.loads(line) for line in (out / "predictions.jsonl").read_text().splitlines()]
    assert len(records) == len(_expected_rows(data))
    manifest = json.loads((out / "run_manifest.json").read_text())
    assert manifest["input_rows"] == len(records)
    known = [r for r in records if r["variant_key"] == "16-2496587-G-C"]
    assert len(known) == 1, "the one valid row must survive its malformed neighbours"
    assert all(m["status"] == "scored_research" for m in known[0]["models"].values())


FILE_LEVEL_REJECTS = {
    "empty": b"",
    "header_only": (HEADER + "\n").encode(),
    "missing_required_column": b"variant_id,chrom,pos,ref,alt\n1,16,2496587,G,C\n",
    "duplicate_required_column": (HEADER + ",variant_id\n" + KNOWN_ROW + ",shadow\n").encode(),
    "not_utf8": (HEADER + "\n" + KNOWN_ROW + "\n").encode() + b"\xff\xfe,GRCh38,1,5,A,G\n",
    "utf16": (HEADER + "\n" + KNOWN_ROW + "\n").encode("utf-16"),
    "field_over_csv_limit": (HEADER + "\n" + KNOWN_ROW + "\nbig,GRCh38,1,5,A," + "G" * 200_000 + "\n").encode(),
    "binary_garbage": bytes(random.Random(3).randrange(256) for _ in range(4096)),
}


@pytest.mark.parametrize("name", sorted(FILE_LEVEL_REJECTS))
def test_cli_file_level_malformations_fail_cleanly_without_partial_output(name, tmp_path, monkeypatch, bundle):
    # UnicodeDecodeError is a ValueError; csv.Error covers the field-size limit.
    with pytest.raises((ValueError, csv.Error)):
        _run(predict_variants, FILE_LEVEL_REJECTS[name], tmp_path, monkeypatch, bundle)
    assert not (tmp_path / "out").exists(), "a rejected input must not leave a partial output directory"


def test_novel_cli_also_rejects_duplicate_required_column(tmp_path, monkeypatch, bundle):
    pytest.importorskip("requests")
    import predict_novel

    with pytest.raises(ValueError, match="Duplicate required column"):
        _run(predict_novel, FILE_LEVEL_REJECTS["duplicate_required_column"], tmp_path, monkeypatch, bundle)
    assert not (tmp_path / "out").exists()
