# Engineering case study: a score needs a defensible contract

## Problem

Random variant splits can allow the same gene on both sides. A high AUC then
answers a different question from generalization to unseen genes. The project
uses connected gene groups, inner-only model/threshold selection and held-out
predictions; the headline AUC is 0.930 versus 0.996 under random splitting.
These differences do not identify a single causal mechanism.

## Statistical repair

Saved predictions allowed an audit without retraining. Whole-group bootstrap
replaced row resampling, the empirical point estimate was separated from the
bootstrap mean, and paired fold differences became descriptive comparisons.
The coding non-truncating subset is not strict missense. The revised intervals
condition on fixed predictions and omit refitting and full CV dependence.
[Reproduction and hashes](statistical_revision.md) document exactly what changed.

## Inference contract

The frozen CLI verifies bundle hashes, excludes held-out gene groups and keeps
failed rows with null scores. The live entrypoint initially accepted a different
conservation score and could confuse GraphQL errors with absent alleles. The
revision blocks incompatible models, validates variant identity, preserves the
frozen absent-frequency transform and treats unknown membership as unknown.
Returning fewer scores is a deliberate consequence of the available evidence.

## How to demonstrate it

1. Run the README frozen demo and inspect variant 48 across five model configurations.
2. Inspect an unavailable input: a null score and reason are an expected result.
3. Run `python -m pytest tests -q` after installing the project, inference
   requirements, `pytest` and `requests`.
4. Follow an empirical estimate and its conditional interval back to the revised
   tables and immutable input hashes.

No new model training, clinical probability calibration or independent external
validation is implied by these engineering checks. Live-service spot checks
remain optional; CI exercises offline contracts and the frozen demonstration.
