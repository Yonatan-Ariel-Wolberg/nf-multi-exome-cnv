# Training module (`--workflow train`)

## Methodology
This module trains an XGBoost classifier from feature matrices and truth labels, producing a serialized model and training report.

## Processes in this module
1. **TRAIN_XGBOOST** – read `*_features.tsv`, join labels, train model, and write outputs.

## Required parameters
- `workflow`: `train`
- `outdir`
- `features_dir`
- `truth_labels`

## Truth-label TSV requirements
The `truth_labels` TSV must include:
- `sample_id`
- `chrom`
- `start`
- `end`
- `cnv_type`
- `truth_label`

Rows are matched by exact keys (`sample_id`, `chrom`, `start`, `end`, `cnv_type`).

## Params template
- `params/params-train.json`
