# InDelible module (`--workflow indelible`)

## Methodology
InDelible identifies structural insertion/deletion events from split-read evidence, builds cohort-level context, annotates candidate events, and supports trio/duo de novo inference.

## Processes in this module
1. **RUN_FETCH** – extract split/clipped read evidence.
2. **RUN_AGGREGATE** – aggregate reads into candidate loci.
3. **RUN_SCORE** – score candidate loci.
4. **RUN_DATABASE** – build cohort-level frequency database.
5. **RUN_ANNOTATE** – annotate with genes/exons and cohort metrics.
6. **RUN_DENOVO_TRIO** – de novo calling for full trios.
7. **RUN_DENOVO_MOM** – de novo calling when only mother is present.
8. **RUN_DENOVO_DAD** – de novo calling when only father is present.
9. **FILTER_INDELIBLE** – filter low-confidence events.
10. **CONVERT_INDELIBLE_TO_VCF** – convert filtered TSV to VCF.
11. **BGZIP_SORT_INDEX_VCF** – sort, compress, and index VCF.
12. **NORMALISE_CNV_QUALITY_SCORES** – normalise caller quality scores.

## Required parameters
- `workflow`: `indelible`
- `outdir`
- `crams`
- `ref`
- `priors`
- `indelible_conf`

## Commonly used optional parameters
- `fai`
- `truth_bed`, `probes_bed` (only when running automatic evaluation in `main.nf`)

## Params template
- `params/params-indelible.json`
