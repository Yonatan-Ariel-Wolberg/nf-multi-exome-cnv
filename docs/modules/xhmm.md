# XHMM module (`--workflow xhmm`)

## Methodology
XHMM calls CNVs from exome read depth using cohort-level depth matrices, PCA-based normalisation, and an HMM for CNV discovery and genotyping.

## Processes in this module
1. **GROUP_BAMS** – partition BAM list into depth-calculation batches.
2. **GATK_DOC** – run GATK DepthOfCoverage per batch.
3. **COMBINE_DOC** – combine depth outputs into one cohort matrix.
4. **CALC_GC_XHMM** – annotate targets with GC information.
5. **FILTER_SAMPLES** – filter outlier targets/samples and center matrix.
6. **RUN_PCA** – run PCA on centered depth matrix.
7. **NORMALISE_PCA** – remove technical variation using PCA components.
8. **FILTER_ZSCORE** – z-score and filter high-variance targets/samples.
9. **FILTER_RD** – align filtered original depth with normalized set.
10. **DISCOVER_CNVS** – discover CNVs.
11. **GENOTYPE_CNVS** – genotype CNVs across samples.
12. **SPLIT_VCF** – split multi-sample VCF into per-sample VCFs.
13. **FILTER_XHMM_CNVS** – apply quality thresholds.
14. **BGZIP_SORT_INDEX_VCF** – sort, compress, and index VCF.
15. **NORMALISE_CNV_QUALITY_SCORES** – normalise caller quality scores.

## Required parameters
- `workflow`: `xhmm`
- `outdir`
- `samplesheet_bams`
- `ref`
- `probes`
- `xhmm_conf`

## Commonly used optional parameters
- `xhmm_batch_size` (default in template: `50`)
- `truth_bed`, `probes_bed` (only when running automatic evaluation in `main.nf`)

## Params template
- `params/params-xhmm.json`
