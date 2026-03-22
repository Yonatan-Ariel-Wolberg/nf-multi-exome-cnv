# CNVkit module (`--workflow cnvkit`)

## Methodology
CNVkit performs target/antitarget coverage analysis, pooled-reference normalisation, and segmentation/calling to produce CNV calls.

## Processes in this module
1. **GENERATE_ACCESS** – derive accessible genome regions.
2. **AUTOBIN** – build target/antitarget bins.
3. **COVERAGE** – compute target and antitarget coverage per sample.
4. **CREATE_POOLED_REFERENCE** – create pooled CNVkit reference.
5. **CALL_CNV** – run fix/segment/call workflow.
6. **EXPORT_RESULTS** – export VCF and BED outputs.
7. **BGZIP_SORT_INDEX_VCF** – sort, compress, and index VCF.
8. **NORMALISE_CNV_QUALITY_SCORES** – normalise caller quality scores.

## Required parameters
- `workflow`: `cnvkit`
- `outdir`
- `bams`
- `fasta`
- `targets`
- `refflat`

## Commonly used optional parameters
- `test_size`
- `test_list`
- `truth_bed`, `probes_bed` (only when running automatic evaluation in `main.nf`)

## Params template
- `params/params-cnvkit.json`
