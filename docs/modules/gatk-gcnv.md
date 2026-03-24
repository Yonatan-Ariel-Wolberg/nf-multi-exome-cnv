# GATK-gCNV module (`--workflow gcnv`)

## Methodology
GATK-gCNV performs cohort-mode germline CNV calling using interval annotations, ploidy modelling, and probabilistic CNV inference with interval sharding.

## Processes in this module
1. **GENERATE_PLOIDY_PRIORS** – generate ploidy priors for contigs.
2. **PREPROCESS_INTERVALS** – preprocess and pad target intervals.
3. **ANNOTATE_INTERVALS** – annotate intervals (e.g., GC content).
4. **COLLECT_READ_COUNTS** – collect read counts into HDF5.
5. **FILTER_INTERVALS** – remove unsuitable intervals.
6. **DETERMINE_PLOIDY_COHORT** – infer ploidy across cohort.
7. **SCATTER_INTERVALS** – split intervals into shards.
8. **GERMLINE_CNV_CALLER_COHORT** – run CNV caller per shard.
9. **POSTPROCESS_CALLS** – combine shard outputs and emit per-sample calls.
10. **BGZIP_SORT_INDEX_VCF** – sort, compress, and index VCF.
11. **NORMALISE_CNV_QUALITY_SCORES** – normalise caller quality scores.

## Required parameters
- `workflow`: `gcnv`
- `outdir`
- `samples_path`
- `fasta`
- `fai`
- `dict`
- `exome_targets`

## Commonly used optional parameters
- `bin_length`
- `padding`
- `is_wgs`
- `scatter_count`
- `truth_bed`, `probes_bed` (only when running automatic evaluation in `main.nf`)

## Params template
- `params/general/params-gatk-gcnv.json`
