# CANOES module (`--workflow canoes`)

## Methodology
CANOES is an exome CNV caller based on read-depth over target probes. In this pipeline it uses cohort-wide read counts, GC-aware normalisation, and the CANOES statistical model to call CNVs before converting output to VCF.

## Processes in this module
1. **CALC_GC_CANOES** – compute GC content per target interval.
2. **BATCH_BAMS** – split BAM inputs into `canoes_batch_size` chunks for parallel depth counting.
3. **GEN_READ_COUNTS** – run `bedtools multicov` per batch/chromosome.
4. **MERGE_READ_COUNTS** – merge per-batch read-count matrices.
5. **RUN_CANOES** – run CANOES CNV calling on merged counts.
6. **FILTER_CANOES_CNVS** – filter low-confidence calls.
7. **CONVERT_CANOES_TO_VCF** – convert CANOES output to VCF.
8. **BGZIP_SORT_INDEX_VCF** – sort, compress, and index VCF.
9. **NORMALISE_CNV_QUALITY_SCORES** – normalise caller quality scores.

## Required parameters
- `workflow`: `canoes`
- `outdir`
- `samplesheet_bams`
- `ref`
- `fai`
- `probes`

## Commonly used optional parameters
- `canoes_batch_size` (default in template: `100`)
- `truth_bed`, `probes_bed` (only when running automatic evaluation in `main.nf`)

## Params template
- `params/params-canoes.json`
