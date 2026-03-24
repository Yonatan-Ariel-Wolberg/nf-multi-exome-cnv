# CLAMMS module (`--workflow clamms`)

## Methodology
CLAMMS uses depth-of-coverage and cohort-derived nearest-neighbour normalisation to train per-sample models and call CNVs.

## Processes in this module
1. **GENERATE_WINDOWS** – create genomic windows with mappability/GC annotations.
2. **SAMTOOLS_DOC** – compute per-sample depth of coverage.
3. **NORMALIZE_DOC** – normalise depth files.
4. **CREATE_PCA_DATA** – run PCA for QC/reference modelling.
5. **GET_PICARD_QC_METRICS** – collect Picard HS metrics.
6. **GET_PICARD_MEAN_INSERT_SIZE** – collect insert-size metrics.
7. **COMBINE_PICARD_QC_METRICS** – combine QC metrics.
8. **CREATE_CUSTOM_REF_PANEL** – build CLAMMS reference panel.
9. **TRAIN_MODELS** – train CLAMMS model per sample.
10. **CALL_CNVS** – call CNVs from trained models.
11. **FILTER_CLAMMS_CNVS** – filter low-confidence calls.
12. **CONVERT_CLAMMS_TO_VCF** – convert calls to VCF.
13. **BGZIP_SORT_INDEX_VCF** – sort, compress, and index VCF.
14. **NORMALISE_CNV_QUALITY_SCORES** – normalise caller quality scores.

## Required parameters
- `workflow`: `clamms`
- `outdir`
- `samplesheet_bams`
- `ref`
- `fai`
- `probes`
- `interval_list`
- `mappability`
- `special_reg`
- `sexinfo`

## Sexinfo file requirements

The `sexinfo` file must be a headerless, tab-separated file with exactly two columns:

- Column 1: `sample_id` – used to look up sex for each sample in the samplesheet.
- Column 2: `sex` – biological sex of the sample; must be `M` (male) or `F` (female).

Example:

```
SAMPLE001	M
SAMPLE002	F
SAMPLE003	M
```

Every sample in the cohort must have an entry in this file. CLAMMS will abort with
an error if a sample ID is not found. The file does not need to contain exactly
the same set of IDs as `samplesheet.tsv` (extra IDs are ignored), but it must
include all sample IDs present in the samplesheet.

## Commonly used optional parameters
- `truth_bed`, `probes_bed` (only when running automatic evaluation in `main.nf`)

## Params template
- `params/general/params-clamms.json`
