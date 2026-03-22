# Feature extraction module (`--workflow feature_extraction`)

## Methodology
This module transforms merged consensus CNV VCFs into per-sample feature matrices for downstream machine learning. It combines consensus metadata with optional caller-specific normalised VCFs and optional sequence/depth annotations.

## Processes in this module
1. **EXTRACT_FEATURES** – generate `{sample}_features.tsv` with structural, support, quality, and optional annotation features from merged CNV calls.

## Required parameters
- `workflow`: `feature_extraction`
- `outdir`
- `merged_vcf_dir`

## Commonly used optional parameters
- `merger_mode` (`survivor` or `truvari`)
- Caller-normalised VCF directories:
  - `canoes_norm_dir`, `clamms_norm_dir`, `xhmm_norm_dir`, `cnvkit_norm_dir`, `gcnv_norm_dir`, `dragen_norm_dir`, `indelible_norm_dir`
- Annotation inputs:
  - `bam_file`
  - `reference_fasta`
  - `bed_file`
  - `mappability_file`
  - `indelible_counts`

## Params templates
- `params/params-survivor-with-features.json` (for SURVIVOR-derived merged VCFs)
- `params/params-truvari-with-features.json` (for Truvari-derived merged VCFs)
