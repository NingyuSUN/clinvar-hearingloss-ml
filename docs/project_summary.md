# Project Summary

This project studies hearing-loss-related variant pathogenicity prediction using ClinVar-derived features and biologically informed machine learning.

The main research problem is whether a model can generalize to variants in genes that were not seen during training. This is important because ordinary random train/test splits can overestimate performance when variants from the same genes appear in both training and testing data.

## Core contribution

The project moves from a simple variant-level classifier to a stricter biological generalization framework:

1. Diagnose random-split overestimation.
2. Evaluate with gene-based train/test separation.
3. Add evolutionary conservation as a cross-gene functional feature.
4. Analyze high-recall operating points using R90.
5. Interpret model behavior with SHAP.

## Key result

Adding conservation improved strict gene-based AUC from 0.653 to 0.812 and reduced R90 false positives from 664 to 394, while maintaining recall at approximately 0.90.

## Current interpretation

SHAP analyses suggest that conservation acts as the primary decision axis when the signal is strong. In an intermediate gray zone, structural features such as deletion, stop-gain, splice-site signal, and variant length become more influential.
