# Evaluation module (`--workflow evaluate`)

## Methodology
This module evaluates CNV calls against a truth BED set by converting per-sample VCFs to BED, combining them, and computing caller performance metrics.

## Processes in this module
1. **VCF_TO_BED** – convert per-sample VCF to BED.
2. **COMBINE_BEDS** – combine all BED files into one cohort BED.
3. **EVALUATE_CALLER** – compute precision/sensitivity/F1 against truth.

## Required parameters
- `workflow`: `evaluate`
- `outdir`
- `vcf_dir`
- `caller`
- `truth_bed`
- `probes_bed`

## Params template
- `params/params-evaluate.json`
