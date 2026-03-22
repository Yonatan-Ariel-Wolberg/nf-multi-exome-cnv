# Normalise module (`--workflow normalise`)

## Methodology
Standalone quality normalisation for pre-existing caller VCF files. This enables post-processing and harmonised quality scaling without rerunning original callers.

## Processes in this module
1. **NORMALISE_CNV_QUALITY_SCORES** – apply caller-specific quality transformation to VCF QUAL scores.

## Required parameters
- `workflow`: `normalise`
- `outdir`
- `vcf_dir`
- `caller` (`CANOES`, `CLAMMS`, `XHMM`, `GCNV`, `CNVKIT`, `DRAGEN`, or `INDELIBLE`)

## Params template
- `params/params-normalise.json`
