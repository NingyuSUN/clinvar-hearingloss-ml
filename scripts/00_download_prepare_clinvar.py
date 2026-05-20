#!/usr/bin/env python
"""Download ClinVar data and prepare hearing-loss variant tables.

This script reproduces the data-origin step of the project:
1. download ClinVar variant_summary.txt.gz from NCBI FTP,
2. filter hearing-loss-related variants using MedGen ID C0018784 or phenotype text,
3. map ClinVar clinical significance into binary labels,
4. save the cleaned table and the early baseline feature matrix.
"""
from __future__ import annotations

import argparse
from pathlib import Path

import numpy as np

from hlpath.clinvar import (
    build_initial_feature_matrix,
    download_variant_summary,
    filter_hearing_loss_variants,
    add_binary_labels,
    keep_grch38,
)


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser()
    parser.add_argument("--raw-dir", default="data/raw", help="Directory for downloaded/raw files.")
    parser.add_argument("--processed-dir", default="data/processed", help="Directory for processed outputs.")
    parser.add_argument("--skip-download", action="store_true", help="Use an existing data/raw/variant_summary.txt.gz file.")
    return parser.parse_args()


def main() -> None:
    args = parse_args()
    raw_dir = Path(args.raw_dir)
    processed_dir = Path(args.processed_dir)
    raw_dir.mkdir(parents=True, exist_ok=True)
    processed_dir.mkdir(parents=True, exist_ok=True)

    variant_summary = raw_dir / "variant_summary.txt.gz"
    if not args.skip_download or not variant_summary.exists():
        print(f"Downloading ClinVar variant summary to {variant_summary} ...")
        download_variant_summary(variant_summary)

    print("Filtering hearing-loss-related ClinVar variants ...")
    hearing = filter_hearing_loss_variants(variant_summary)
    hearing_path = raw_dir / "hearing_loss_clinvar_clean.csv"
    hearing.to_csv(hearing_path, index=False)

    print("Creating binary labels from ClinicalSignificance ...")
    clean = add_binary_labels(hearing)
    clean_path = processed_dir / "hearing_loss_df_clean.csv"
    clean.to_csv(clean_path, index=False)

    print("Saving GRCh38-only table for external annotation ...")
    grch38 = keep_grch38(clean)
    grch38_path = processed_dir / "hearing_loss_grch38.csv"
    grch38.to_csv(grch38_path, index=False)

    print("Building early baseline feature matrix ...")
    X, y, _ = build_initial_feature_matrix(clean)
    X.to_parquet(processed_dir / "hearing_loss_X.parquet")
    np.save(processed_dir / "hearing_loss_y.npy", y)

    print("Done.")
    print(f"Raw filtered table:       {hearing_path}")
    print(f"Clean labeled table:      {clean_path}")
    print(f"GRCh38 annotation input:  {grch38_path}")
    print(f"Baseline feature matrix:  {processed_dir / 'hearing_loss_X.parquet'}")


if __name__ == "__main__":
    main()
