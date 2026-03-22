# nf-multicaller-exomecnv

A Nextflow pipeline for comprehensive genomic analysis combining multiple variant callers with exome and copy number variation (CNV) detection. This workflow is based on a workflow by **Dr Phelelani T Mpangase** ([phelelani](https://github.com/phelelani)): [nf-exomecnv](https://github.com/phelelani/nf-exomecnv), which incorporated CANOES, CLAMMS, XHMM and InDelible, but stopped at the actual variant calling component, without including structural variant merging and *in silico* validation with a ML classifier. The ```modules-icav2-dragen.nf``` module is based off a script by **Regan Cannell** ([RWCannell](https://github.com/RWCannell)): [cram_input_dragen_ica_workflow](https://github.com/SBIMB/ica-elwazi/tree/main/nextflow_workflows/cram_input_dragen_ica_workflow). The ```modules-cnvkit.nf``` module follows the instructions laid out in the [CNVkit Copy Number Calling Pipeline](https://cnvkit.readthedocs.io/en/stable/pipeline.html). The ```modules-gatk-gcnv.nf``` module follows the instructions laid out in the [GATK Guide on How to Call Rare Germline Copy Number Variants](https://gatk.broadinstitute.org/hc/en-us/articles/360035531152--How-to-Call-rare-germline-copy-number-variants).

## Tools

**CANOES**: [GitHub](https://github.com/ShenLab/CANOES); [DOI](https://doi.org/10.1093/nar/gku345)

**CLAMMS**: [GitHub](https://github.com/rgcgithub/clamms); [DOI](https://doi.org/10.1093/bioinformatics/btv547)

**XHMM**: [GitHub](https://github.com/RRafiee/XHMM); [DOI](https://doi.org/10.1002/0471142905.hg0723s81)

**InDelible**: [GitHub](https://github.com/HurlesGroupSanger/indelible); [DOI](https://doi.org/10.1016/j.ajhg.2021.09.010)

**ICAv2**: [GitHub](https://github.com/umccr/illumination/tree/v2); [ICAv2 Commands](https://help.ica.illumina.com/command-line-interface/cli-indexcommands)

**Illumina's DRAGEN Germline Enrichment**: [Description of DRAGEN Enrichment](https://www.illumina.com/products/by-type/informatics-products/basespace-sequence-hub/apps/dragen-enrichment.html); [DOI](https://doi.org/10.1038/s41587-024-02382-1)

**The Broad Institute's GATK-gCNV**: [GitHub](https://github.com/broadinstitute/gatk); [DOI](https://doi.org/10.1038/s41588-023-01449-0)

**CNVkit**: [GitHub](https://github.com/etal/cnvkit); [DOI](https://doi.org/10.1371/journal.pcbi.1004873)

**Truvari**: [GitHub](https://github.com/ACEnglish/truvari); [DOI](https://doi.org/10.1186/s13059-022-02840-6)

**SURVIVOR**: [GitHub](https://github.com/fritzsedlazeck/SURVIVOR); [DOI](https://doi.org/10.1038/ncomms14061)

## Overview

This pipeline implements a robust bioinformatics workflow for analyzing genomic data using:
- Multiple variant calling approaches
- Exome analysis
- Copy number variation (CNV) detection

Built with Nextflow for reproducible, scalable, and portable data analysis.

## Workflow Steps

The pipeline is structured as a sequence of modular stages.  Each stage can be
run independently; in a typical end-to-end analysis the stages are executed in
the order described below.

### Step 1 — Individual CNV callers

Seven independent CNV callers are provided.  Run one or more of them first to
produce per-sample, per-caller VCF files that feed into the downstream
consensus and ML stages.

#### CANOES (`--workflow canoes`)
1. **CALC_GC_CANOES** – Compute GC content per target region using GATK.
2. **BATCH_BAMS** – Split the cohort BAM list into fixed-size batches
   (`canoes_batch_size`) to bound memory and file-descriptor use.
3. **GEN_READ_COUNTS** – Run `bedtools multicov` per batch and chromosome to
   produce per-batch read-count matrices.
4. **MERGE_READ_COUNTS** – Merge per-batch matrices into a single cohort-wide
   read-count matrix.
5. **RUN_CANOES** – Normalise read counts by mean coverage, select and weight
   reference samples using non-negative least squares (NNLS), and estimate
   per-probe variance with a GAM (GC content as covariate).  CNVs are then
   called using a **Hidden Markov Model (HMM) in which the Negative Binomial
   distribution serves as the emission model**: for each probe, `dnbinom` log-
   probabilities are computed under the deletion (CN=1), normal (CN=2), and
   duplication (CN=3) states and passed directly into the Viterbi decoder — the
   Negative Binomial and the HMM are therefore a single integrated model, not
   separate steps.  The Phred-scaled quality score Q_SOME (probability that the
   CNV event exists) is computed via the forward–backward algorithm using the
   same emission probabilities.
6. **FILTER_CANOES_CNVS** – Apply quality thresholds to remove low-confidence
   calls.
7. **CONVERT_CANOES_TO_VCF** – Convert the filtered CSV output to VCF format,
   annotating each record with `TOOL=CANOES`.

#### XHMM (`--workflow xhmm`)
1. **GROUP_BAMS** – Partition the BAM list into groups of `xhmm_batch_size` for
   parallel GATK `DepthOfCoverage` jobs.
2. **GATK_DOC** – Run GATK `DepthOfCoverage` per BAM group to produce
   per-group sample-interval depth summaries.
3. **COMBINE_DOC** – Merge all per-group depth-of-coverage files into a single
   cohort-wide read-depth matrix using `xhmm --mergeGATKdepths`.
4. **CALC_GC_XHMM** – Annotate each target interval with its GC content using
   GATK `AnnotateIntervals` and identify extreme-GC targets (GC < 10 % or
   > 90 %) for exclusion.
5. **FILTER_SAMPLES** – Remove targets that fail size or coverage thresholds and
   samples that fail mean read-depth thresholds, then mean-centre the retained
   per-target read depths using `xhmm --matrix --centerData`.
6. **RUN_PCA** – Perform principal component analysis on the mean-centred depth
   matrix using `xhmm --PCA` to capture systematic technical variation.
7. **NORMALISE_PCA** – Remove the top PCA components (PVE_mean factor 0.7) from
   the mean-centred data to produce a PCA-normalised read-depth matrix using
   `xhmm --normalize`.
8. **FILTER_ZSCORE** – Re-filter the PCA-normalised matrix to remove
   high-variance targets (maxSdTargetRD > 30) and compute per-sample z-scores
   using `xhmm --matrix --zScoreData`, emitting lists of excluded targets and
   samples.
9. **FILTER_RD** – Intersect the original read-depth matrix with the combined
   target and sample exclusion lists produced by steps 5 and 8 so that the
   original and normalised matrices cover identical loci and samples.
10. **DISCOVER_CNVS** – Discover CNVs across the cohort by running XHMM's HMM
    on the z-score matrix (`xhmm --discover`), using the original filtered
    read depths as the reference.
11. **GENOTYPE_CNVS** – Genotype all discovered CNVs in every sample using
    `xhmm --genotype` and emit a combined multi-sample VCF.
12. **SPLIT_VCF** – Split the multi-sample VCF into individual per-sample VCFs,
    retaining only non-reference, non-missing genotype calls, and add the `chr`
    prefix to chromosome names if absent.
13. **FILTER_XHMM_CNVS** – Apply quality-score filters (EQ, SQ, and NDQ
    thresholds) to each per-sample VCF using `bcftools filter` to remove
    low-confidence CNV calls.

#### CLAMMS (`--workflow clamms`)
1. **GENERATE_WINDOWS** – Define fixed-size genomic windows over the target
   intervals.
2. **SAMTOOLS_DOC** – Compute per-sample depth-of-coverage with `samtools`.
3. **NORMALIZE_DOC** – Apply CLAMMS nearest-neighbour reference-panel
   normalisation to each sample's coverage.
4. **CREATE_PCA_DATA** – Run PCA on normalised coverages for sample QC.
5. **GET_PICARD_QC_METRICS** / **GET_PICARD_MEAN_INSERT_SIZE** – Collect Picard
   HS metrics and mean insert size per sample.
6. **COMBINE_PICARD_QC_METRICS** – Aggregate per-sample Picard metrics into a
   cohort-wide table.
7. **CREATE_CUSTOM_REF_PANEL** – Build the cohort reference panel from all
   normalised coverage files.
8. **TRAIN_MODELS** – Train a per-sample CLAMMS model against the reference
   panel.
9. **CALL_CNVS** – Call CNVs per sample using the trained model.
10. **FILTER_CLAMMS_CNVS** – Apply quality filters.
11. **CONVERT_CLAMMS_TO_VCF** – Convert filtered BED output to VCF, annotating
    each record with `TOOL=CLAMMS`.

#### GATK-gCNV (`--workflow gcnv`)
1. **GENERATE_PLOIDY_PRIORS** – Prepare cohort-level ploidy priors.
2. **PREPROCESS_INTERVALS** – Standardise and bin the target intervals.
3. **ANNOTATE_INTERVALS** – Attach GC-content and mappability annotations to
   each interval.
4. **COLLECT_READ_COUNTS** – Count reads per interval per sample.
5. **FILTER_INTERVALS** – Remove low-quality or low-coverage intervals.
6. **DETERMINE_PLOIDY_COHORT** – Estimate per-sample contig ploidy across the
   cohort.
7. **SCATTER_INTERVALS** – Partition intervals into shards for parallel calling
   (`scatter_count`).
8. **GERMLINE_CNV_CALLER_COHORT** – Run the GATK probabilistic cohort-mode CNV
   model per shard.
9. **POSTPROCESS_CALLS** – Merge shards, apply posterior filters, and produce
   per-sample VCFs annotated with `TOOL=GCNV`.

#### CNVkit (`--workflow cnvkit`)
1. **GENERATE_ACCESS** – Identify accessible (mappable) regions of the genome.
2. **AUTOBIN** – Optimise target and antitarget bin sizes for the cohort.
3. **COVERAGE** – Compute per-sample target and antitarget coverage.
4. **CREATE_POOLED_REFERENCE** – Build a pooled reference from all sample
   coverages for GC and bias correction.
5. **CALL_CNV** – Call copy number per sample against the pooled reference.
6. **EXPORT_RESULTS** – Export calls to VCF and BED, annotating VCF records
   with `TOOL=CNVKIT`.

#### DRAGEN (`--workflow dragen`)
1. **UPLOAD_CRAM_FILES** – Upload CRAM and CRAI files for each sample to the
   ICAv2 cloud platform and record the resulting ICAv2 file IDs.
2. **GET_STATIC_FILES** – Append the ICAv2 IDs of the reference genome and
   annotation files to the combined upload manifest.
3. **CHECK_FILE_STATUS** – Poll ICAv2 to confirm that every uploaded file is
   available before submitting the analysis job.
4. **START_ANALYSIS_BATCH** – Submit a DRAGEN Germline Enrichment batch job to
   ICAv2 using the verified file manifest.
5. **CHECK_ANALYSIS_STATUS** – Poll ICAv2 until the DRAGEN analysis job reaches
   a terminal state (succeeded or failed).
6. **DOWNLOAD_ANALYSIS_OUTPUT** – Retrieve the result VCFs, JSON reports, and
   BAMs from ICAv2.
7. **DELETE_DATA** – Remove temporary CRAM files and the output folder from
   ICAv2 to avoid storage costs.
8. **ADD_DRAGEN_TOOL_ANNOTATION** – Annotate each downloaded VCF record with
   `TOOL=DRAGEN` in the INFO field.

#### InDelible (`--workflow indelible`)
1. **RUN_FETCH** – Extract split and soft-clipped reads from each CRAM/BAM file
   using `indelible.py fetch` to collect per-read evidence for indel events.
2. **RUN_AGGREGATE** – Aggregate per-read split-read information to a
   position-level view of the data using `indelible.py aggregate`, combining
   BAM alignment context and clipped-read evidence.
3. **RUN_SCORE** – Score each candidate position based on read information and
   local sequence context using `indelible.py score`.
4. **RUN_DATABASE** – Construct a cohort-wide allele-frequency and breakpoint
   database from all scored files using `indelible.py database`, providing the
   population-level context required for annotation.
5. **RUN_ANNOTATE** – Enrich each scored file with gene and exon annotations
   and merge in allele-frequency and breakpoint database results using
   `indelible.py annotate`.
6. **RUN_DENOVO_TRIO** – Identify de novo indel mutations in complete trios
   (proband, mother, father) using `indelible.py denovo` with both parental
   BAMs.
7. **RUN_DENOVO_MOM** – Identify de novo mutations when only the mother's BAM
   is available.
8. **RUN_DENOVO_DAD** – Identify de novo mutations when only the father's BAM
   is available.
9. **FILTER_INDELIBLE** – Remove low-confidence calls by applying
   population-frequency thresholds on the allele-frequency (`AF_freq`) and
   breakpoint-frequency (`BP_freq`) columns of the annotated TSV.
10. **CONVERT_INDELIBLE_TO_VCF** – Convert the filtered TSV to VCF format using
    `indelible_tsv_to_vcf.py`, annotating each record with `TOOL=INDELIBLE`.

---

### Step 2 — VCF post-processing (all callers)

After each caller produces its VCF, two shared processes are applied
automatically:

1. **BGZIP_SORT_INDEX_VCF** – Block-compress (`bgzip`), coordinate-sort
   (`bcftools sort`), and index (`tabix`) the VCF, producing a sorted
   `.vcf.gz` + `.tbi` pair.
2. **NORMALISE_CNV_QUALITY_SCORES** – Rescale the caller-native quality score
   to a standardised 0–100 scale and write it to the `QUAL` field, storing
   the original score in the `OQ` FORMAT field.

---

### Step 3 — Consensus calling

Once two or more callers have been run, their VCF directories are passed to one
of the two consensus modules.  Both produce a per-sample consensus call set by
merging the independent caller outputs.

#### SURVIVOR (`--workflow survivor`)
1. **RUN_SURVIVOR_MERGE (union)** – Merge all per-caller VCFs with
   `min_support = 1` (any caller) and a merge distance of 1,000 bp to produce
   a **union** call set containing every variant seen by at least one caller.
2. **RUN_SURVIVOR_MERGE (intersection)** – Re-merge with `min_support = 2`
   (two or more callers) to produce an **intersection** call set retaining only
   variants supported by ≥ 2 callers.

#### Truvari (`--workflow truvari`)
1. **MERGE_VCFS** – Concatenate and sort all per-caller VCFs into a single
   multi-sample VCF.
2. **COLLAPSE_VCFS** – Run `truvari collapse` with size-reciprocal-overlap
   (`pctsize = 0.5`) and breakpoint-overlap (`pctovl = 0.5`) thresholds to
   cluster redundant calls and produce a single collapsed consensus VCF.

---

### Step 4 — Feature extraction (`--workflow feature_extraction`)

The feature-extraction stage transforms the merged consensus VCF into a
structured feature matrix suitable for machine learning.

1. **EXTRACT_FEATURES** (`bin/feature_extraction.py`) – For each variant in
   the merged VCF, compute ~40 features drawn from multiple input sources
   (see the table below).

The output is a per-sample TSV (`{sample}_features.tsv`) for use in training
and prediction.

#### Feature table

The table below lists every column written to the feature matrix.
The **Used in training** column reflects which columns are retained by
`bin/train_xgboost.py` (non-numeric and pure identifier columns
`sample_id`, `chrom`, and `cnv_type` are dropped before model fitting;
`start` and `end` are numeric and therefore kept).

**Structural / location features**

| Feature | Description | Source file / input | Used in training |
|---|---|---|---|
| `chrom` | Chromosome name (string, e.g. `chr1`) | Merged SURVIVOR/Truvari VCF | No (identifier) |
| `start` | CNV start position (0-based) | Merged SURVIVOR/Truvari VCF | Yes |
| `end` | CNV end position | Merged SURVIVOR/Truvari VCF | Yes |
| `chrom_encoded` | Chromosome as integer: autosomes 1–22, X→23, Y→24, other→0 | Derived from `chrom` | Yes |
| `cnv_size` | CNV length in base-pairs | Derived from `start`/`end` | Yes |
| `size_label` | Ordinal size-bin (1–12) matching CN-Learn's `SIZE_LABEL` scheme | Derived from `cnv_size` | Yes |
| `cnv_type` | CNV type indicator: 1 = DUP, 0 = DEL/other | Merged SURVIVOR/Truvari VCF | No (non-numeric string in some formats) |

**Concordance**

| Feature | Description | Source file / input | Used in training |
|---|---|---|---|
| `concordance` | Number of callers supporting this event | Merged SURVIVOR/Truvari VCF (`SUPP` INFO field or `MatchId` counts) | Yes |

**Per-caller support flags**

| Feature | Description | Source file / input | Used in training |
|---|---|---|---|
| `is_canoes` | 1 if CANOES called this event, else 0 | Merged VCF `SUPP_VEC` / `TOOL` INFO field | Yes |
| `is_clamms` | 1 if CLAMMS called this event, else 0 | Merged VCF `SUPP_VEC` / `TOOL` INFO field | Yes |
| `is_xhmm` | 1 if XHMM called this event, else 0 | Merged VCF `SUPP_VEC` / `TOOL` INFO field | Yes |
| `is_gatk_gcnv` | 1 if GATK-gCNV called this event, else 0 | Merged VCF `SUPP_VEC` / `TOOL` INFO field | Yes |
| `is_cnvkit` | 1 if CNVkit called this event, else 0 | Merged VCF `SUPP_VEC` / `TOOL` INFO field | Yes |
| `is_dragen` | 1 if DRAGEN called this event, else 0 | Merged VCF `SUPP_VEC` / `TOOL` INFO field | Yes |
| `is_indelible` | 1 if INDELIBLE called this event, else 0 | Merged VCF `SUPP_VEC` / `TOOL` INFO field | Yes |

**Per-caller normalised quality scores**

Raw caller-native scores are converted to a common 0–100 scale by
`bin/normalise_cnv_caller_quality_scores.py` before feature extraction.
The conversion formula for each caller is:

| Caller | Raw score | Conversion |
|---|---|---|
| CANOES | `Q_SOME` (FORMAT) | `min(100, 100 × Q_SOME / 80)` |
| CLAMMS | `Q_SOME` (FORMAT) | `min(100, 100 × Q_SOME / 500)` |
| XHMM | `SQ` (FORMAT) | `min(100, 100 × SQ / 60)` |
| GATK-gCNV | `QS` (FORMAT) | `min(100, 100 × QS / threshold)` |
| CNVkit | `CNQ` (FORMAT) | `min(100, 100 × CNQ / 20)` |
| DRAGEN | native `QUAL` | linear rescaling to 0–100 |
| INDELIBLE | synthetic Phred | `SR_TOTAL × (AVG_MAPQ / 60) × 100` |

| Feature | Description | Source file / input | Used in training |
|---|---|---|---|
| `qual_norm_canoes` | QUAL_norm score for CANOES (0–100) | Normalised CANOES VCF (`QUAL` field) | Yes |
| `qual_norm_clamms` | QUAL_norm score for CLAMMS (0–100) | Normalised CLAMMS VCF (`QUAL` field) | Yes |
| `qual_norm_xhmm` | QUAL_norm score for XHMM (0–100) | Normalised XHMM VCF (`QUAL` field) | Yes |
| `qual_norm_gatk_gcnv` | QUAL_norm score for GATK-gCNV (0–100) | Normalised GATK-gCNV VCF (`QUAL` field) | Yes |
| `qual_norm_cnvkit` | QUAL_norm score for CNVkit (0–100) | Normalised CNVkit VCF (`QUAL` field) | Yes |
| `qual_norm_dragen` | QUAL_norm score for DRAGEN (0–100) | Normalised DRAGEN VCF (`QUAL` field) | Yes |
| `qual_norm_indelible` | QUAL_norm score for INDELIBLE (0–100) | Normalised INDELIBLE VCF (`QUAL` field) | Yes |

**Aggregate quality summaries**

| Feature | Description | Source file / input | Used in training |
|---|---|---|---|
| `max_qual_norm` | Maximum QUAL_norm across all callers supporting this event | Derived from per-caller normalised VCFs | Yes |
| `mean_qual_norm_supported` | Mean QUAL_norm across callers that support this event (NaN when concordance = 0) | Derived from per-caller normalised VCFs | Yes |

**Genomic annotations**

| Feature | Description | Source file / input | Used in training |
|---|---|---|---|
| `n_probes` | Number of capture-target probes overlapping the CNV (NaN if no BED supplied) | Capture-target BED file | Yes |
| `n_probes_flank` | Number of capture-target probes in the flanking region (NaN if no BED supplied) | Capture-target BED file | Yes |
| `rd_ratio` | Mean read depth inside the CNV divided by mean read depth in flanking regions (NaN if no BAM/CRAM supplied) | BAM/CRAM alignment file | Yes |
| `gc_content` | GC nucleotide fraction of the CNV interval (NaN if no FASTA supplied) | Reference FASTA | Yes |
| `mappability` | Weighted-mean mappability score over the CNV interval (NaN if no mappability file supplied) | Mappability BED file | Yes |

**Probe-level log2-ratio statistics**

Computed per capture probe from BAM read depths normalised against the flanking-probe baseline.
All three features are NaN when the capture BED or BAM/CRAM is absent.

| Feature | Description | Source file / input | Used in training |
|---|---|---|---|
| `l2r_mean` | Mean log2(probe depth / flanking baseline) across CNV probes; negative for DEL, positive for DUP | BAM/CRAM + capture-target BED | Yes |
| `l2r_dev` | Variance (not standard deviation) of per-probe log2 ratios within the CNV; high values indicate noisy or false calls | BAM/CRAM + capture-target BED | Yes |
| `l2r_flank_diff` | `l2r_mean` minus the mean flanking-probe log2 ratio; captures deviation from the immediate neighbourhood | BAM/CRAM + capture-target BED | Yes |

**Caller-specific secondary INFO metrics**

| Feature | Description | Source file / input | Used in training |
|---|---|---|---|
| `xhmm_rd` | XHMM read-depth Z-score (NaN if XHMM not run) | Normalised XHMM VCF (`RD` INFO field) | Yes |
| `cnvkit_weight` | CNVkit bin weight (NaN if CNVkit not run) | Normalised CNVkit VCF (`WEIGHT` INFO field) | Yes |
| `cnvkit_log2` | CNVkit log2 copy-ratio (NaN if CNVkit not run) | Normalised CNVkit VCF (`LOG2` INFO field) | Yes |
| `dragen_sm` | DRAGEN sample median (NaN if DRAGEN not run) | Normalised DRAGEN VCF (`SM` INFO field) | Yes |
| `dragen_sd` | DRAGEN sample standard deviation (NaN if DRAGEN not run) | Normalised DRAGEN VCF (`SD` INFO field) | Yes |

**INDELIBLE split-read metrics**

| Feature | Description | Source file / input | Used in training |
|---|---|---|---|
| `total_sr` | Total split-read count supporting the event (NaN if INDELIBLE not run) | INDELIBLE counts TSV | Yes |
| `sr_entropy` | Shannon entropy of split-read positions (NaN if INDELIBLE not run) | INDELIBLE counts TSV | Yes |
| `mapq_avg` | Average mapping quality of split reads (NaN if INDELIBLE not run) | INDELIBLE counts TSV | Yes |
| `dual_split` | 1 if both breakends have split-read support, else 0 (NaN if INDELIBLE not run) | INDELIBLE counts TSV | Yes |

---

### Step 5 — ML training (`--workflow train`)

1. **TRAIN_XGBOOST** (`bin/train_xgboost.py`) – Load all
   `*_features.tsv` files alongside a truth-label TSV, balance the training
   set using SMOTE oversampling, and train an XGBoost gradient-boosted
   classifier with stratified k-fold cross-validation.  Outputs:
   - `cnv_model.json` – the serialised XGBoost model.
   - `training_report.txt` – cross-validation precision, recall, and F1-score.

**Truth-label TSV requirements (including truthset CNVs labelled as true)**

The `--truth_labels` file must be a tab-delimited TSV containing at least:

- `sample_id`
- `chrom`
- `start`
- `end`
- `cnv_type` (e.g., `DEL` / `DUP`)
- `truth_label`

Rows are joined to the feature matrix using **exact** matches on
`sample_id`, `chrom`, `start`, `end`, and `cnv_type` (inner join).

For truthset CNVs that should be treated as **true positives** during
training, set `truth_label = 1`. False positives should use
`truth_label = 0`.

Example:

```tsv
sample_id	chrom	start	end	cnv_type	truth_label
SAMPLE_001	chr1	100000	130000	DEL	1
SAMPLE_001	chr1	500000	520000	DUP	0
```

---

### Step 6 — Performance evaluation (`--workflow evaluate`)

1. **VCF_TO_BED** (`bin/vcf_to_bed.py`) – Convert each per-sample VCF to a
   5-column BED file (CHR, START, STOP, CNV_TYPE, SAMPLE_ID).
2. **COMBINE_BEDS** – Concatenate all per-sample BED files into a single
   cohort-wide call-set BED.
3. **EVALUATE_CALLER** (`bin/evaluate_caller_performance.py`) – Compare the
   unified call set against a truth BED at probe level and report precision,
   sensitivity, and F1-score per caller.

---

### Step 7 — Standalone normalisation (`--workflow normalise`)

For cohorts where per-caller VCFs already exist (e.g. produced outside this
pipeline), the normalisation stage can be run in isolation:

1. **NORMALISE_CNV_QUALITY_SCORES** – Accept a directory of raw caller VCFs
   and a caller name (CANOES, CLAMMS, XHMM, GCNV, CNVKIT, DRAGEN, or
   INDELIBLE), and produce bgzipped, sorted, indexed VCFs with standardised
   quality scores.

---

## Features

- **Multi-caller approach**: Integrates multiple variant calling algorithms for robust variant detection
- **Exome analysis**: Specialized processing for whole exome sequencing (WES) data
- **CNV detection**: Comprehensive copy number variation analysis
- **Nextflow-based**: Ensures reproducibility and scalability across different computing environments
- **Containerized**: Includes Docker/Singularity support for dependency management

## Requirements

### Operating System
- **Linux** or **Unix** (macOS, HPC clusters running Linux)
  > Windows is **not** supported. Use WSL2 or a remote Linux server if you are on Windows.

### Software
- **Apptainer** (formerly Singularity): >= 1.0  
  All tool containers are pulled and executed via Apptainer. Install it following the [official guide](https://apptainer.org/docs/admin/main/installation.html). ([DOI](https://doi.org/10.1371/journal.pone.0177459))
- **Nextflow**: >= 20.04  ([DOI](https://doi.org/10.1038/nbt.3820)) 
- Install with:
  ```bash
  curl -s https://get.nextflow.io | bash
  ```
- **Java**: >= 8 (required by Nextflow)  
  Install OpenJDK:
  - Debian/Ubuntu: `sudo apt install default-jdk`
  - RHEL/CentOS/Rocky: `sudo yum install java-1.8.0-openjdk`
  - macOS (Homebrew): `brew install openjdk`

### Libraries and Dependencies
All bioinformatics tools (CANOES, CLAMMS, XHMM, GATK, CNVkit, InDelible, ICAv2 CLI, SURVIVOR, Truvari, bcftools, bedtools, Picard) are provided as pre-built container images and are downloaded automatically by Apptainer at runtime. No manual tool installation is required beyond Apptainer and Nextflow themselves.

### Computing Resources
- Adequate CPU cores (parallel execution recommended)
- Sufficient disk space for intermediate and output files
- Memory requirements depend on dataset size

## Installation

1. Clone this repository:
```bash
git clone https://github.com/Yonatan-Ariel-Wolberg/nf-multicaller-exomecnv.git
cd nf-multicaller-exomecnv
```

2. Build the icav2_cli_v2.43.0.sif
```bash
apptainer build icav2_cli_v2.43.0.sif icav2.def
```

## Usage

Each workflow is run by providing `--workflow <name>` along with a params file. Edit the corresponding `params/*.json` file to point to your actual input files before running.

### Individual CNV callers

The individual callers should be run first. Their VCF outputs feed into the consensus modules (SURVIVOR / Truvari).

```bash
# CANOES — exome CNV caller using read-depth, NNLS reference weighting and Negative Binomial HMM
nextflow run main.nf --workflow canoes -params-file params/params-canoes.json

# XHMM — exome CNV caller using GATK depth-of-coverage and PCA
nextflow run main.nf --workflow xhmm -params-file params/params-xhmm.json

# CLAMMS — exome CNV caller using nearest-neighbour normalisation
nextflow run main.nf --workflow clamms -params-file params/params-clamms.json

# GATK-gCNV — cohort-mode germline CNV caller from the Broad Institute
nextflow run main.nf --workflow gcnv -params-file params/params-gatk-gcnv.json

# CNVkit — coverage-based CNV caller with pooled reference normalisation
nextflow run main.nf --workflow cnvkit -params-file params/params-cnvkit.json

# DRAGEN — Illumina DRAGEN germline enrichment pipeline via ICAv2
nextflow run main.nf --workflow dragen -params-file params/params-icav2-dragen.json

# InDelible — small insertion/deletion caller from split reads
nextflow run main.nf --workflow indelible -params-file params/params-indelible.json
```

### Consensus modules (run after individual callers)

Point the `*_dir` paths in the params file to the VCF output directories produced by the callers above, then run:

```bash
# SURVIVOR — merge per-caller VCFs into a consensus call set
nextflow run main.nf --workflow survivor -params-file params/params-survivor.json

# Truvari — collapse per-caller VCFs into a consensus call set
nextflow run main.nf --workflow truvari -params-file params/params-truvari.json
```

### Parameter files

All required and optional parameters for each workflow are documented in the corresponding `params/*.json` template:

| Params file | Workflow |
|---|---|
| `params/params-canoes.json` | `--workflow canoes` |
| `params/params-xhmm.json` | `--workflow xhmm` |
| `params/params-clamms.json` | `--workflow clamms` |
| `params/params-gatk-gcnv.json` | `--workflow gcnv` |
| `params/params-cnvkit.json` | `--workflow cnvkit` |
| `params/params-icav2-dragen.json` | `--workflow dragen` |
| `params/params-indelible.json` | `--workflow indelible` |
| `params/params-survivor.json` | `--workflow survivor` |
| `params/params-truvari.json` | `--workflow truvari` |

## Running on the Wits UI Cluster

The workflow has been tested on the **University of the Witwatersrand (Wits) UI HPC cluster**, which runs **Linux** and uses **SLURM** as its job scheduler.  Use the `wits` profile to activate the cluster-specific settings.

### Quick start on the Wits cluster

```bash
# Load the required modules (adjust module names to match the cluster's module system)
module load nextflow
module load apptainer   # or: module load singularity (apptainer is already loaded by default on the Wits UI Cluster)

# Run a workflow — replace <workflow> and the params file as needed
nextflow run main.nf \
    -profile wits,medium \
    --workflow canoes \
    -params-file params/params-canoes.json
```

### What the `wits` profile does

| Setting | Value |
|---|---|
| `executor` | `slurm` |
| `process.queue` | `batch` |
| `process.cpus` | 4 |
| `process.memory` | 16 GB |
| `process.time` | 72 h |
| `executor.queueSize` | 500 |
| `executor.submitRateLimit` | 10/sec |

### Cluster-specific bind mounts

The global `singularity.runOptions` in `nextflow.config` binds several Wits file-system paths into every container:

```
-B /dataB/aux
-B /dataG/ddd
-B /dataG/ddd-2023
-B /home/ywolberg
```

Edit `params.runOptions` in `nextflow.config` (or pass `--runOptions` on the command line) if your data lives elsewhere on the cluster.

### GATK / Picard constraint

GATK and Picard processes request nodes with AVX2 support (`--constraint=avx2`) and exclude node `n04` (`--exclude=n04`).  Adjust `process.clusterOptions` in `nextflow.config` if the Wits cluster topology changes.

### ICAv2 / DRAGEN on the Wits cluster

The DRAGEN workflow submits runs to the ICAv2 cloud platform from a local container.  Before running `--workflow dragen`:

1. Build the ICAv2 CLI container (if not already done):
   ```bash
   apptainer build icav2_cli_v2.43.0.sif icav2.def
   ```
2. Update the container path in `nextflow.config` (`withLabel: 'icav2-dragen'`) to the location where you stored the `.sif` file.
3. Ensure your ICA credentials are present at `~/.icav2/` on the cluster head node so they can be bind-mounted into the container.

## Scalability: running 50, 300, or 2000 samples

The pipeline ships with three cohort-size profiles that pre-configure the
scheduler, concurrency, and per-caller batch sizes for small, medium, and large
runs.  Combine a cohort-size profile with the site profile using `-profile`:

```bash
# ≤ 50 samples
nextflow run main.nf -profile wits,small  -params-file params/params-canoes.json

# 50–300 samples (default behaviour if no cohort-size profile is specified)
nextflow run main.nf -profile wits,medium -params-file params/params-canoes.json

# 300–2000 samples
nextflow run main.nf -profile wits,large  -params-file params/params-canoes.json
```

### What each profile controls

| Setting | `small` (≤50) | `medium` (50–300) | `large` (300–2000) |
|---|---|---|---|
| `executor.queueSize` | 100 | 500 | 2000 |
| `executor.submitRateLimit` | 5/sec | 10/sec | 20/sec |
| `process.maxForks` | 25 | 100 | 200 |
| `canoes_batch_size` | 50 | 100 | 200 |
| `xhmm_batch_size` | 25 | 50 | 100 |
| GATK memory | 16 GB | 16 GB | 32 GB |
| `large_mem` label | 64 GB / 8 CPUs | 64 GB / 8 CPUs | 128 GB / 16 CPUs |
| `gpu_or_high_cpu` label | 32 GB / 16 CPUs | 32 GB / 16 CPUs | 64 GB / 32 CPUs |

### Caller-specific notes

#### CANOES (`--workflow canoes`)
CANOES uses `bedtools multicov` to count reads across all BAMs simultaneously.
Running a single `multicov` job over thousands of BAMs at once exhausts both
memory and file-descriptor limits.  The pipeline therefore splits the BAM list
into fixed-size batches (`canoes_batch_size`) and runs one `multicov` job per
batch per chromosome; the per-batch matrices are then merged into a single
cohort-wide matrix **before** the CANOES R script normalises read counts,
selects reference samples via NNLS, and calls CNVs using an HMM in which the
Negative Binomial distribution is the integrated emission model (the NB
log-probabilities are fed directly into the Viterbi decoder — they are a
single combined model, not separate steps).  The batch size therefore controls
**only parallelism and
memory use**, not normalisation accuracy — every sample contributes to the
final normalisation regardless of batch size.

- **`canoes_batch_size`** (default `100`; override in params file or via profile):
  - Up to 50 samples → set to `50` (1 batch; minimal overhead).
  - Up to 300 samples → set to `100` (3 batches per chromosome).
  - Up to 2000 samples → set to `200` (10 batches per chromosome).

#### XHMM (`--workflow xhmm`)
XHMM uses GATK `DepthOfCoverage`, which can require significant memory per job
when given a large BAM list.  The BAM list is split into fixed-size groups
(`xhmm_batch_size`) and one `GATK_DOC` job is submitted per group.  All groups
are merged by `xhmm --mergeGATKdepths` before PCA normalisation, so the batch
size controls **only job memory**, not the normalisation cohort.

- **`xhmm_batch_size`** (default `50`; override in params file or via profile):
  - Up to 50 samples → set to `25`–`50` (1–2 batches).
  - Up to 300 samples → set to `50` (up to 6 batches).
  - Up to 2000 samples → set to `100` (up to 20 batches; keeps memory per job manageable).

#### CLAMMS (`--workflow clamms`)
All CLAMMS steps (depth-of-coverage, normalisation, model training, CNV calling)
are fully per-sample and run in parallel automatically via Nextflow's data-flow
model.  No batch-size tuning is needed; `process.maxForks` limits the number of
concurrent per-sample jobs.  The custom reference panel is built from **all**
normalised coverage files in the cohort, ensuring accurate nearest-neighbour
normalisation regardless of cohort size.

#### GATK-gCNV (`--workflow gcnv`)
gCNV trains a cohort-level model and is naturally designed for large cohorts.
Interval parallelism is controlled by `scatter_count` in
`params/params-gatk-gcnv.json` (default `5000`).  For the cohort-level steps
(`FilterIntervals`, `DetermineGermlineContigPloidy`, `GermlineCNVCaller`) that
collect all samples at once, increase the memory/CPU limits using the `large`
profile when running hundreds or thousands of samples.

#### CNVkit (`--workflow cnvkit`)
CNVkit creates a pooled reference from **all** sample coverages before calling
CNVs, so every sample contributes to bias correction and GC normalisation.
Concurrency is governed by `process.maxForks`.  Use `test_size` (integer) or
`test_list` (comma-separated sample IDs) in the params file to restrict a run
to a subset of samples during development.

#### DRAGEN / ICAv2 (`--workflow dragen`)
DRAGEN runs are submitted to the ICAv2 platform.  `maxUploadForks` (default `4`)
limits simultaneous CRAM uploads.  Increase this if your network bandwidth and
ICA account limits allow it.

#### SURVIVOR / Truvari (`--workflow survivor` / `--workflow truvari`)
These consensus modules are fast and lightweight regardless of cohort size; no
additional tuning is required.

### Minimum cohort sizes for accurate normalisation

Each CNV caller's normalisation method requires a minimum number of samples to
produce statistically reliable results.  Running with fewer samples than the
recommended minimum will not cause the pipeline to fail, but CNV call quality
will be reduced.

| Caller | Normalisation method | Recommended minimum | Notes |
|--------|---------------------|--------------------:|-------|
| **CANOES** | Mean-coverage normalisation + NNLS reference weighting + Negative Binomial HMM | 30 | NNLS reference-panel quality degrades with very small cohorts; 30+ samples recommended for reliable reference weighting |
| **XHMM** | Mean-centering + PCA | 30 | XHMM's own filters may exclude all targets/samples with very small cohorts; 50+ strongly recommended |
| **CLAMMS** | Nearest-neighbour reference panel | 30 | Nearest-neighbour selection quality degrades when the reference pool is small |
| **GATK-gCNV** | Probabilistic cohort model | 30 | The Broad recommends ≥ 30 samples for the COHORT mode model to be well-calibrated |
| **CNVkit** | Pooled reference (bias correction) | 10 | A pooled reference with as few as 10 samples provides reasonable GC and mappability correction; larger cohorts improve accuracy |
| **DRAGEN** | Illumina platform normalisation | N/A | Normalisation is handled internally by the DRAGEN pipeline on the ICAv2 platform. ≥ 5 samples (inclusive) are required to create the 'in-run' PoN. ≥ 30 samples required for Targeted CNV Calling |

> **Note:** `canoes_batch_size` and `xhmm_batch_size` control how many BAMs are
> processed per parallel job, not how many samples participate in normalisation.
> All per-batch results are merged into a single cohort-wide matrix before any
> normalisation step runs, so these parameters have no effect on normalisation
> accuracy.

### Disk and I/O considerations
- Intermediate files (depth-of-coverage, read-count matrices, normalised
  coverages) grow proportionally with cohort size.  Ensure the working directory
  has sufficient free space: budget ~2–5 GB per sample.
- Setting `process.stageInMode = 'symlink'` (the default) avoids copying large
  BAM/CRAM files into each work directory.
