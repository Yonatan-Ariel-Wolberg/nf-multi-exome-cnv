# Common VCF post-processing module (`modules/common/modules-common.nf`)

## Methodology
Most caller modules run two shared post-processing steps: VCF standardisation (sort/compress/index and tool annotation) and quality-score normalisation.

## Processes in this module
1. **BGZIP_SORT_INDEX_VCF** – ensure caller VCFs are bgzipped, coordinate-sorted, tabix-indexed, and tagged with `TOOL` metadata.
2. **NORMALISE_CNV_QUALITY_SCORES** – convert caller-specific quality scores to a common scale.

## Required parameters (when reused by caller modules)
- `outdir`
- caller label suffix/path inputs passed from upstream module outputs

## Notes
This is an internal shared module used by caller workflows (CANOES, XHMM, CLAMMS, CNVkit, GATK-gCNV, DRAGEN, InDelible).
