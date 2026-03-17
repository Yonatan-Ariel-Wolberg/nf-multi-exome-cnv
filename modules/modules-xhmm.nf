#!/usr/bin/env nextflow
nextflow.enable.dsl=2

// =====================================================================================
// XHMM MODULE
// =====================================================================================

// =====================================================================================
// GLOBAL FILE INSTANTIATION
// =====================================================================================
params.xhmm_batch_size = 50  // BAMs per GATK DepthOfCoverage job; increase for fewer but larger jobs
ref       = file(params.ref, type: 'file')
probes    = file(params.probes, type: 'file')
xhmm_conf = file(params.xhmm_conf, type: 'file')
outdir    = file(params.outdir, type: 'dir')

// =====================================================================================
// SPLIT BAM FILES TO GROUPS
// =====================================================================================
process GROUP_BAMS {
    tag { 'group_bams' }

    input:
    path(bams)

    output:
    path("bam_group_*"), emit: bam_groups
    
    script:
    def batch_size = params.xhmm_batch_size as int
    """
    sort ${bams} > bam_list_sorted.txt
    split -l ${batch_size} bam_list_sorted.txt --numeric-suffixes=000 --suffix-length=3 --additional-suffix=.list bam_group_
    """
}

// RUN GATK FOR DEPTH OF COVERAGE (FOR SAMPLES IN EACH GROUP):
process GATK_DOC {
    tag { "${group}" }
    label 'gatk'
    publishDir "${outdir}/out_XHMM", mode: 'copy', overwrite: true
    
    input:
    tuple val(group), path(list)
    
    output:
    tuple val(group), path("${group}.DATA.sample_interval_summary"), emit: bam_group_doc

    script:
    mem = task.memory.toGiga() - 4   // Reserve memory for GATK processing
    
    """
    gatk --java-options "-XX:+UseSerialGC -Xss456k -Xms2g -Xmx${mem}g -Djava.io.tmpdir=\$PWD" \
        DepthOfCoverage \
        -I ${list} \
        -L ${probes} \
        -R ${ref} \
        --max-depth-per-sample 5000 \
        --verbosity INFO \
        --omit-depth-output-at-each-base true \
        --omit-locus-table true \
        --min-base-quality 0 \
        --read-filter MappingQualityReadFilter \
        --minimum-mapping-quality 20 \
        --start 1 --stop 5000 --nBins 200 \
        --include-ref-n-sites true \
        --count-type COUNT_READS \
        --output ${group}.DATA
    """
}

// COMBINES GATK DEPTH-OF-COVERAGE OUTPUTS FOR MULTIPLE SAMPLES (AT SAME LOCI):
process COMBINE_DOC {
    tag { 'combine_doc' }
    label 'xhmm'
    publishDir "${outdir}/out_XHMM", mode: 'copy', overwrite: true
    
    input:
    path(list)
    
    output:
    path("DATA.RD.txt"), emit: combined_doc
    
    script:
    """
    for i in *.sample_interval_summary; do sed 's/,/	/g' \$i > \${i%.sample_interval_summary}.fixed_sample_interval_summary ; done
    ls *.fixed_sample_interval_summary > the_list
    xhmm --mergeGATKdepths -o DATA.RD.txt --GATKdepthsList the_list
    """
}

// OPTIONALLY, RUN GATK TO CALCULATE THE PER-TARGET GC CONTENT AND CREATE A LIST OF EXTREME GC CONTENT TARGETS:
process CALC_GC_XHMM {
    tag { 'calc_gc' }
    label 'gatk'
    publishDir "${outdir}/out_XHMM", mode: 'copy', overwrite: true
    
    output:
    path("DATA.locus_GC.txt"), emit: gc_targets
    path("extreme_gc_targets.txt"), emit: extreme_gc_targets
    
    script:
    """
    gatk AnnotateIntervals \
        -R ${ref} \
        -L ${probes} \
        --interval-merging-rule OVERLAPPING_ONLY \
        -O gc.tmp

    grep -v "^@" gc.tmp | grep -v "^CONTIG" | awk '{ print \$1":"\$2"-"\$3"\t"\$4 }' > DATA.locus_GC.txt
    cat DATA.locus_GC.txt | awk '{ if (\$2 < 0.1 || \$2 > 0.9) print \$1 }' > extreme_gc_targets.txt
    """
}

// FILTERS SAMPLES AND TARGETS AND THEN MEAN-CENTERS THE TARGETS:
process FILTER_SAMPLES {
    tag { 'filter_samples' }
    label 'xhmm'
    publishDir "${outdir}/out_XHMM", mode: 'copy', overwrite: true    

    input:
    path(combined_doc)
    path(extreme_gc_targets)
    
    output:
    path("DATA.filtered_centered.RD.txt"), emit: filtered_centered
    path("DATA.filtered_centered.RD.txt.filtered_targets.txt"), emit: excluded_filtered_targets
    path("DATA.filtered_centered.RD.txt.filtered_samples.txt"), emit: excluded_filtered_samples
    
    script:
    """
    xhmm --matrix -r DATA.RD.txt --centerData --centerType target \
        -o DATA.filtered_centered.RD.txt \
        --outputExcludedTargets DATA.filtered_centered.RD.txt.filtered_targets.txt \
        --outputExcludedSamples DATA.filtered_centered.RD.txt.filtered_samples.txt \
        --excludeTargets extreme_gc_targets.txt \
        --minTargetSize 10 --maxTargetSize 10000 \
        --minMeanTargetRD 10 --maxMeanTargetRD 500 \
        --minMeanSampleRD 25 --maxMeanSampleRD 200 \
        --maxSdSampleRD 150
   """
}

// RUNS PCA ON MEAN-CENTERED DATA:
process RUN_PCA {
    tag { 'run_pca' }
    label 'xhmm'
    publishDir "${outdir}/out_XHMM", mode: 'copy', overwrite: true

    input:
    path(filtered_centered)
    
    output:
    tuple val("pca_data"), path("DATA.RD_PCA*"), emit: pca_data

    script:
    """
    xhmm --PCA -r DATA.filtered_centered.RD.txt --PCAfiles DATA.RD_PCA
    """
}

// NORMALIZES MEAN-CENTERED DATA USING PCA INFORMATION:
process NORMALISE_PCA {
    tag { 'norm_pca' }
    label 'xhmm'
    publishDir "${outdir}/out_XHMM", mode: 'copy', overwrite: true

    input:
    path(filtered_centered)
    tuple val(pca), path(pca_data)
    
    output:
    path("DATA.PCA_normalized.txt"), emit: data_pca_norm

    script:
    """
    xhmm --normalize -r DATA.filtered_centered.RD.txt \
        --PCAfiles DATA.RD_PCA \
        --normalizeOutput DATA.PCA_normalized.txt \
        --PCnormalizeMethod PVE_mean --PVE_mean_factor 0.7
    """
}

// FILTERS AND Z-SCORE CENTERS (BY SAMPLE) THE PCA-NORMALIZED DATA:
process FILTER_ZSCORE {
    tag { 'filter_zscore' }
    label 'xhmm'
    publishDir "${outdir}/out_XHMM", mode: 'copy', overwrite: true
    
    input:
    path(data_pca_norm)
    
    output:
    path("DATA.PCA_normalized.filtered.sample_zscores.RD.txt"), emit: pca_norm_zscore
    path("DATA.PCA_normalized.filtered.sample_zscores.RD.txt.filtered_targets.txt"), emit: excluded_zscore_targets
    path("DATA.PCA_normalized.filtered.sample_zscores.RD.txt.filtered_samples.txt"), emit: excluded_zscore_samples
    
    script:
    """
    xhmm --matrix -r DATA.PCA_normalized.txt \
        --centerData --centerType sample --zScoreData \
        -o DATA.PCA_normalized.filtered.sample_zscores.RD.txt \
        --outputExcludedTargets DATA.PCA_normalized.filtered.sample_zscores.RD.txt.filtered_targets.txt \
        --outputExcludedSamples DATA.PCA_normalized.filtered.sample_zscores.RD.txt.filtered_samples.txt \
        --maxSdTargetRD 30
    """
}

// FILTERS ORIGINAL READ-DEPTH DATA TO BE THE SAME AS FILTERED, NORMALIZED DATA:
process FILTER_RD {
    tag { 'filter_rd' }
    label 'xhmm'
    publishDir "${outdir}/out_XHMM", mode: 'copy', overwrite: true

    input:
    path(combined_doc)
    path(excluded_filtered_targets)
    path(excluded_zscore_targets)
    path(excluded_filtered_samples)
    path(excluded_zscore_samples)
    
    output:
    path("DATA.same_filtered.RD.txt"), emit: orig_filtered

    script:
    """
    xhmm --matrix -r DATA.RD.txt \
        --excludeTargets DATA.filtered_centered.RD.txt.filtered_targets.txt \
        --excludeTargets DATA.PCA_normalized.filtered.sample_zscores.RD.txt.filtered_targets.txt \
        --excludeSamples DATA.filtered_centered.RD.txt.filtered_samples.txt \
        --excludeSamples DATA.PCA_normalized.filtered.sample_zscores.RD.txt.filtered_samples.txt \
        -o DATA.same_filtered.RD.txt
    """
}

// DISCOVERS CNVS IN NORMALIZED DATA:
process DISCOVER_CNVS {
    tag { 'discover_cnvs' }
    label 'xhmm'
    publishDir "${outdir}/out_XHMM", mode: 'copy', overwrite: true

    input:
    path(orig_filtered)
    path(pca_norm_zscore)
    
    output:
    path("*"), emit: cnvs
    
    script:
    """
    xhmm --discover -p ${xhmm_conf} \
        -r DATA.PCA_normalized.filtered.sample_zscores.RD.txt \
        -R DATA.same_filtered.RD.txt \
        -c DATA.xcnv -a DATA.aux_xcnv -s DATA
    """
}

// GENOTYPES DISCOVERED CNVS IN ALL SAMPLES:
process GENOTYPE_CNVS {
    tag { 'genotype_cnvs' }
    label 'xhmm'
    publishDir "${outdir}/out_XHMM", mode: 'copy', overwrite: true
    
    input:
    path(orig_filtered)
    path(pca_norm_zscore)
    path(cnvs)
    
    output:
    path("*"), emit: genotypes
    
    script:
    """
    xhmm --genotype -p ${xhmm_conf} \
        -r DATA.PCA_normalized.filtered.sample_zscores.RD.txt \
        -R DATA.same_filtered.RD.txt \
        -g DATA.xcnv -F ${ref} \
        -v DATA.vcf
    """
}

// SPLITS COMBINED VCF INTO INDIVIDUAL SAMPLE VCFs
process SPLIT_VCF {
    tag { 'split_vcf' }
    label 'bcftools'
    publishDir "${outdir}/out_XHMM", mode: 'copy', overwrite: true
    
    input:
    path(vcf_file)
    
    output:
    path("*_XHMM_output.vcf"), emit: individual_vcfs
    
    script:
    """
    # Split the multi-sample VCF by sample, keeping only variants actually
    # called for each sample (GT != ref and GT != missing), and ensure the
    # CHROM column carries the 'chr' prefix required by downstream tools.
    for sample in \$(bcftools query -l ${vcf_file}); do
        bcftools view -s \${sample} ${vcf_file} | \
        bcftools view -i 'GT!="." && GT!="0"' | \
        awk 'BEGIN{OFS="\t"} /^#/{print; next} {if (\$1 !~ /^chr/) \$1="chr"\$1; print}' \
        > \${sample}_XHMM_output.vcf
    done
    """
}

// FILTERS CNV OUTPUTS USING BCFTOOLS BASED ON EQ, SQ, AND NDQ VALUES
process FILTER_XHMM_CNVS {
    tag { 'filter_cnvs' }
    label 'bcftools'
    publishDir "${outdir}/out_XHMM", mode: 'copy', overwrite: true
    
    input:
    path(individual_vcfs)
    
    output:
    path("*_filtered.vcf"), emit: filtered_cnvs
    
    script:
    """
    # Loop over the individual VCF files and apply bcftools filter
    for vcf in *.vcf; do
        bcftools filter -e 'FORMAT/EQ<=60 || FORMAT/SQ<=60 || FORMAT/NDQ<=60' \${vcf} -o \${vcf%.vcf}_filtered.vcf
    done
    """
}

// Process to compress, sort, index, and annotate each XHMM VCF with TOOL=XHMM
process BGZIP_SORT_INDEX_VCF {
    tag "${vcf_file.simpleName}"
    label 'bcftools'
    publishDir "${outdir}/out_XHMM/vcfs", mode: 'copy', overwrite: true

    input:
    path vcf_file

    output:
    path("*.sorted.vcf.gz"),     emit: sorted_vcf
    path("*.sorted.vcf.gz.tbi"), emit: sorted_vcf_index

    script:
    def sample_name = vcf_file.simpleName
    def sorted_gz   = "${sample_name}.sorted.vcf.gz"
    """
    # Create extra header line for the TOOL INFO field
    printf '##INFO=<ID=TOOL,Number=1,Type=String,Description="Calling tool">\\n' > extra_header.txt

    # Build a BED annotation file with TOOL=XHMM for every variant
    bcftools query -f '%CHROM\\t%POS0\\t%END\\n' ${vcf_file} | \\
        awk 'BEGIN{OFS="\\t"} {print \$1, \$2, \$3, "XHMM"}' | \\
        bgzip -c > ${sample_name}_tool_annot.bed.gz
    tabix -p bed ${sample_name}_tool_annot.bed.gz

    # Annotate the VCF with TOOL=XHMM in the INFO field
    bcftools annotate \\
        -a ${sample_name}_tool_annot.bed.gz \\
        -c CHROM,FROM,TO,INFO/TOOL \\
        -h extra_header.txt \\
        ${vcf_file} \\
        -O v -o ${sample_name}_annotated.vcf

    # Compress the annotated VCF with bgzip
    bgzip -c ${sample_name}_annotated.vcf > ${sample_name}_annotated.vcf.gz

    # Sort the compressed VCF
    bcftools sort ${sample_name}_annotated.vcf.gz -o ${sorted_gz} -O z

    # Index the sorted VCF
    tabix -p vcf ${sorted_gz}

    # Remove intermediate files
    rm -f ${sample_name}_tool_annot.bed.gz ${sample_name}_tool_annot.bed.gz.tbi \\
          ${sample_name}_annotated.vcf ${sample_name}_annotated.vcf.gz
    """
}

// =====================================================================================
// SUB-WORKFLOW TO CHAIN THE PROCESSES TOGETHER
// =====================================================================================

workflow XHMM {
    take:
    bam_list_ch  // path: file containing BAM file paths, one per line

    main:
    // Step 1: Group BAMs into batches for parallel depth-of-coverage calculation
    GROUP_BAMS(bam_list_ch)

    bam_groups_ch = GROUP_BAMS.out.bam_groups
        .flatten()
        .map { f -> tuple(f.simpleName, f) }

    // Step 2: Calculate depth of coverage per BAM group
    GATK_DOC(bam_groups_ch)

    // Step 3: Combine all group depth-of-coverage outputs
    COMBINE_DOC(GATK_DOC.out.bam_group_doc.map { group, file -> file }.collect())

    // Step 4: Calculate per-target GC content and identify extreme GC targets
    CALC_GC_XHMM()

    // Step 5: Filter samples and targets, then mean-center the data
    FILTER_SAMPLES(COMBINE_DOC.out.combined_doc, CALC_GC_XHMM.out.extreme_gc_targets)

    // Step 6: Run PCA on mean-centered data
    RUN_PCA(FILTER_SAMPLES.out.filtered_centered)

    // Step 7: Normalize mean-centered data using PCA information
    NORMALISE_PCA(FILTER_SAMPLES.out.filtered_centered, RUN_PCA.out.pca_data)

    // Step 8: Filter and z-score center the PCA-normalized data
    FILTER_ZSCORE(NORMALISE_PCA.out.data_pca_norm)

    // Step 9: Filter original read-depth data to match the normalized filtered data
    FILTER_RD(
        COMBINE_DOC.out.combined_doc,
        FILTER_SAMPLES.out.excluded_filtered_targets,
        FILTER_ZSCORE.out.excluded_zscore_targets,
        FILTER_SAMPLES.out.excluded_filtered_samples,
        FILTER_ZSCORE.out.excluded_zscore_samples
    )

    // Step 10: Discover CNVs in the normalized data
    DISCOVER_CNVS(FILTER_RD.out.orig_filtered, FILTER_ZSCORE.out.pca_norm_zscore)

    // Step 11: Genotype discovered CNVs across all samples
    GENOTYPE_CNVS(
        FILTER_RD.out.orig_filtered,
        FILTER_ZSCORE.out.pca_norm_zscore,
        DISCOVER_CNVS.out.cnvs
    )

    // Step 12: Split the combined VCF into per-sample VCFs
    vcf_ch = GENOTYPE_CNVS.out.genotypes
        .flatten()
        .filter { it.name == 'DATA.vcf' }

    SPLIT_VCF(vcf_ch)

    // Step 13: Filter per-sample CNV calls by quality scores
    FILTER_XHMM_CNVS(SPLIT_VCF.out.individual_vcfs)

    // Step 14: Compress, sort, index, and annotate each VCF with TOOL=XHMM
    BGZIP_SORT_INDEX_VCF(FILTER_XHMM_CNVS.out.filtered_cnvs.flatten())

    emit:
    sorted_vcf       = BGZIP_SORT_INDEX_VCF.out.sorted_vcf
    sorted_vcf_index = BGZIP_SORT_INDEX_VCF.out.sorted_vcf_index
}
