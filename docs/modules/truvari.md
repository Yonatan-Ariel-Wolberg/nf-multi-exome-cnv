# Truvari module (`--workflow truvari`)

## Methodology
Truvari consensus calling first merges per-caller VCFs per sample, then collapses redundant calls based on size and overlap thresholds.

## Processes in this module
1. **MERGE_VCFS** – combine per-caller VCFs into merged per-sample VCFs.
2. **COLLAPSE_VCFS** – collapse overlapping/redundant variants into consensus calls.

## Required parameters
- `workflow`: `truvari`
- `outdir`
- At least two caller VCF directories from:
  - `canoes_dir`
  - `clamms_dir`
  - `xhmm_dir`
  - `cnvkit_dir`
  - `gcnv_dir`
  - `dragen_dir`
  - `indelible_dir`

## Params template
- `params/params-truvari.json`
