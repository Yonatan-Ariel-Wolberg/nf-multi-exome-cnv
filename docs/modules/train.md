# Training module (`--workflow train`)

## Methodology
This module trains an XGBoost classifier from feature matrices and truth labels, producing a serialized model, training report, ROC/PR curves, and SHAP feature-attribution artifacts.

## Processes in this module
1. **TRAIN_XGBOOST** – read `*_features.tsv`, join labels, train model, and write outputs.

## Container
- `withLabel: 'train'` uses `docker://quay.io/biocontainers/xgboost:2.0.3--py310h4aa3b51_0`.

## Required parameters
- `workflow`: `train`
- `outdir`
- `features_dir`
- `truth_labels`

## Optional parameters
- `probes_bed` – capture-target BED used for fallback matching by shared probes.
- `min_shared_probes` – minimum shared-probe count for fallback matching (default: `1`).

## Truth-label TSV requirements
The `truth_labels` TSV must include:
- `sample_id`
- `chrom`
- `start`
- `end`
- `cnv_type`
- `truth_label`

Rows are first matched by exact keys (`sample_id`, `chrom`, `start`, `end`, `cnv_type`).
When `probes_bed` is provided, any still-unmatched rows are then matched by shared
capture probes on the same `sample_id`, `chrom`, and `cnv_type`.

## Params template
- `params/params-train.json`

## Outputs
- `cnv_model.json` – trained XGBoost model.
- `training_report.txt` – summary of sample counts, class balance, match counts, and AUC metrics.
- `roc_curve.svg` / `roc_curve.tsv` – ROC plot and underlying points.
- `pr_curve.svg` / `pr_curve.tsv` – precision-recall plot and underlying points.
- `shap_values.tsv` – per-call SHAP values for numeric model features.
- `shap_summary_bar.svg` – SHAP mean absolute contribution bar chart.
- `shap_summary_beeswarm.svg` – SHAP beeswarm plot for top features.
