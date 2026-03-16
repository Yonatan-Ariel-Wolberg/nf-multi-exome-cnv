#!/usr/bin/env nextflow
nextflow.enable.dsl=2

// =====================================================================================
// CLAMMS MODULE
// =====================================================================================

// =====================================================================================
// GLOBAL FILE INSTANTIATION
// =====================================================================================
params.ref           = ''
params.probes        = ''
params.interval_list = ''
params.mappability   = ''
params.special_reg   = ''
params.sexinfo       = ''
params.outdir        = 'results'

// Define file variables for inputs
ref           = file(params.ref)
probes        = file(params.probes)
interval_list = file(params.interval_list)
mappability   = file(params.mappability)
special_reg   = file(params.special_reg)
sexinfo       = file(params.sexinfo)
outdir        = file(params.outdir)

// =====================================================================================
// PROCESSES FOR CLAMMS
// =====================================================================================

// Process to generate windows for analysis
process generateWindows {
    tag "Generate_Windows"
    label 'clamms|bedtools'

    output:
    path "windows.bed", emit: windows
    
    script:
    """
    #!/bin/bash
    set -euo pipefail
    echo "Generating windows..."
    annotate_windows.sh ${ref} ${probes} ${mappability} ${special_reg} windows.bed
    echo "Finished generating windows: windows.bed"
    """
}

// Process to calculate depth of coverage using samtools
process samtoolsDOC {
    tag "${sample_id}"
    label 'clamms|bedtools'

    input:
    tuple val(sample_id), path(bam), path(bai)
    path windows

    output:
    tuple val(sample_id), path("${sample_id}.coverage.bed"), emit: coverage
    
    script:
    """
    #!/bin/bash
    set -euo pipefail
    echo "Calculating depth of coverage for ${sample_id}..."
    samtools bedcov ${windows} ${bam} > ${sample_id}.coverage.bed
    echo "Coverage file created for ${sample_id}: ${sample_id}.coverage.bed"
    """
}

// Process to normalize depth of coverage
process normalizeDOC {
    tag "${sample_id}"
    label 'clamms|bedtools'

    input:
    tuple val(sample_id), path(coverage)
    path windows

    output:
    tuple val(sample_id), path("${sample_id}.norm.coverage.bed"), emit: norm_coverage
    
    script:
    """
    #!/bin/bash
    set -euo pipefail
    echo "Normalizing depth of coverage for ${sample_id}..."
    normalize_coverage ${coverage} ${windows} > ${sample_id}.norm.coverage.bed
    echo "Normalized coverage file created for ${sample_id}: ${sample_id}.norm.coverage.bed"
    """
}

// Process to create PCA input data
process createPCAData {
    tag "PCA_Data"
    label 'clamms|bedtools'

    input:
    path norm_covs

    output:
    path "pca_data.txt", emit: pca_data
    
    script:
    """
    #!/bin/bash
    set -euo pipefail
    echo "Creating PCA data..."
    ls *.norm.coverage.bed > norm_covs.list
    echo "PCA matrix placeholder" > pca_data.txt
    echo "PCA data created: pca_data.txt"
    """
}

// Process to get QC metrics using Picard
process getPicardQCMetrics {
    tag "${sample_id}"
    label 'picard'

    input:
    tuple val(sample_id), path(bam), path(bai)

    output:
    tuple val(sample_id), path("${sample_id}.hs_metrics.txt"), emit: qc_metrics
    
    script:
    """
    #!/bin/bash
    set -euo pipefail
    echo "Collecting QC metrics for ${sample_id}..."
    picard CollectHsMetrics \\
        -I ${bam} \\
        -O ${sample_id}.hs_metrics.txt \\
        -R ${ref} \\
        -BI ${interval_list} \\
        -TI ${interval_list}
    echo "QC metrics collected: ${sample_id}.hs_metrics.txt"
    """
}

// Process to get mean insert size metrics using Picard
process getPicardMeanInsertSize {
    tag "${sample_id}"
    label 'picard'

    input:
    tuple val(sample_id), path(bam), path(bai)

    output:
    tuple val(sample_id), path("${sample_id}.insert_size_metrics.txt"), emit: ins_size_metrics
    
    script:
    """
    #!/bin/bash
    set -euo pipefail
    echo "Collecting mean insert size metrics for ${sample_id}..."
    picard CollectInsertSizeMetrics \\
        -I ${bam} \\
        -O ${sample_id}.insert_size_metrics.txt \\
        -H ${sample_id}.insert_size_histogram.pdf
    echo "Insert size metrics collected: ${sample_id}.insert_size_metrics.txt"
    """
}

// Process to combine all QC metrics from Picard
process combinePicardQCMetrics {
    tag "Combine_QC"
    label 'clamms|bedtools'

    input:
    path metrics

    output:
    path "combined_qc_metrics.txt", emit: qcs_metrics
    
    script:
    """
    #!/bin/bash
    set -euo pipefail
    echo "Combining QC metrics..."
    combine_picard_qc_metrics.sh > combined_qc_metrics.txt
    echo "Combined QC metrics saved to: combined_qc_metrics.txt"
    """
}

// Process to create a custom reference panel
process createCustomRefPanel {
    tag "Ref_Panel"
    label 'clamms|bedtools'

    input:
    path norm_covs
    path pca_data
    path combined_qc
    
    output:
    path "*.ref.panel.files", emit: ref_panel
    
    script:
    """
    #!/bin/bash
    set -euo pipefail
    echo "Creating custom reference panel..."
    custom_ref_panel.sh ${combined_qc} ${pca_data}
    echo "Custom reference panel created."
    """
}

// Process to train models using the reference panel
process trainModels {
    tag "${sample_id}"
    label 'clamms|bedtools'

    input:
    tuple val(sample_id), path(ref_panel_file)
    path windows
    path norm_covs

    output:
    tuple val(sample_id), path("${sample_id}.models.bed"), emit: sample_models
    
    script:
    """
    #!/bin/bash
    set -euo pipefail
    echo "Training models for ${sample_id}..."
    fit_models ${windows} ${ref_panel_file} > ${sample_id}.models.bed
    echo "Models trained for ${sample_id}: ${sample_id}.models.bed"
    """
}

// Process to call CNVs from the trained models
process callCNVs {
    tag "${sample_id}"
    label 'clamms|bedtools'
    publishDir "${outdir}/out_CLAMMS/calls", mode: 'copy', overwrite: true

    input:
    tuple val(sample_id), path(sample_models), path(norm_cov)

    output:
    tuple val(sample_id), path("${sample_id}_CLAMMS_calls.bed"), emit: cnvs
    
    script:
    """
    #!/bin/bash
    set -euo pipefail
    echo "Calling CNVs for ${sample_id}..."
    call_cnv ${sample_models} ${norm_cov} > ${sample_id}_CLAMMS_calls.bed
    echo "CNVs called for ${sample_id}: ${sample_id}_CLAMMS_calls.bed"
    """
}

// Process to filter CNVs based on criteria
process filterCLAMMSCNVs {
    tag "Filter_CLAMMS"
    label 'clamms|bedtools'
    publishDir "${outdir}/out_CLAMMS", mode: 'copy', overwrite: true

    input:
    path cnvs

    output:
    path "Sample_CNVs_filtered.bed", emit: filtered_cnvs
    
    script:
    """
    #!/bin/bash
    set -euo pipefail
    echo "Filtering CNVs..."
    cat ${cnvs} > Sample_CNVs.bed
    awk '\$9 >= 500 && \$10 >= 0' Sample_CNVs.bed > Sample_CNVs_filtered.bed || true
    echo "Filtered CNVs saved: Sample_CNVs_filtered.bed"
    """
}


// Process to convert CLAMMS BED to VCF using the Python script
process convertClammsToVcf {
    tag "BED_TO_VCF"
    label 'clamms|bedtools'
    publishDir "${outdir}/out_CLAMMS/vcfs", mode: 'copy', overwrite: true

    input:
    path input_file
    path sample_file
    path fai_file

    output:
    path("*_CLAMMS_output.vcf"), emit: vcfs

    script:
    """
    python3 clamms_bed_to_vcf.py --input_file ${input_file} --sample_file ${sample_file} --output_dir . --fai_file ${fai_file}
    echo "Conversion from BED to VCF completed"
    """
}

// Process to bgzip, sort, index, and annotate the produced VCF file with TOOL=CLAMMS
process BGZIP_SORT_INDEX_VCF {
    tag "${vcf_file.simpleName}"
    label 'bcftools'
    publishDir "${outdir}/out_CLAMMS/vcfs", mode: 'copy', overwrite: true

    input:
    path vcf_file

    output:
    path("*.sorted.vcf.gz"), emit: sorted_vcf
    path("*.sorted.vcf.gz.tbi"), emit: sorted_vcf_index

    script:
    def sample_name = vcf_file.simpleName
    def sorted_gz   = "${sample_name}.sorted.vcf.gz"
    """
    # Create extra header line for the TOOL INFO field
    printf '##INFO=<ID=TOOL,Number=1,Type=String,Description="Calling tool">\\n' > extra_header.txt

    # Build a BED annotation file with TOOL=CLAMMS for every variant
    bcftools query -f '%CHROM\\t%POS0\\t%END\\n' ${vcf_file} | \\
        awk 'BEGIN{OFS="\\t"} {print \$1, \$2, \$3, "CLAMMS"}' | \\
        bgzip -c > ${sample_name}_tool_annot.bed.gz
    tabix -p bed ${sample_name}_tool_annot.bed.gz

    # Annotate the VCF with TOOL=CLAMMS in the INFO field
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

workflow CLAMMS {
    take:
    bam_ch         // channel: tuple(sample_id, bam, bai)
    fai_ch         // path: reference FASTA index (.fai)
    sample_file_ch // path: sample file for VCF conversion

    main:
    // Step 1: Generate analysis windows
    generateWindows()

    // Step 2: Calculate depth of coverage per sample
    samtoolsDOC(bam_ch, generateWindows.out.windows)

    // Step 3: Normalize depth of coverage per sample
    normalizeDOC(samtoolsDOC.out.coverage, generateWindows.out.windows)

    // Step 4: Collect normalized coverages and create PCA data
    createPCAData(normalizeDOC.out.norm_coverage.map { sample_id, norm_cov -> norm_cov }.collect())

    // Step 5: Get Picard QC metrics per sample
    getPicardQCMetrics(bam_ch)
    getPicardMeanInsertSize(bam_ch)

    // Step 6: Combine all Picard QC metrics
    all_metrics_ch = getPicardQCMetrics.out.qc_metrics
        .mix(getPicardMeanInsertSize.out.ins_size_metrics)
        .map { sample_id, metrics -> metrics }
        .collect()

    combinePicardQCMetrics(all_metrics_ch)

    // Step 7: Create custom reference panel using all normalized coverages
    createCustomRefPanel(
        normalizeDOC.out.norm_coverage.map { sample_id, norm_cov -> norm_cov }.collect(),
        createPCAData.out.pca_data,
        combinePicardQCMetrics.out.qcs_metrics
    )

    // Step 8: Train models per sample using its reference panel file
    ref_panel_per_sample_ch = createCustomRefPanel.out.ref_panel
        .flatten()
        .map { file -> tuple(file.name.replace('.ref.panel.files', ''), file) }

    trainModels(
        ref_panel_per_sample_ch,
        generateWindows.out.windows,
        normalizeDOC.out.norm_coverage.map { sample_id, norm_cov -> norm_cov }.collect()
    )

    // Step 9: Call CNVs per sample by joining models with normalized coverage
    call_cnv_input_ch = trainModels.out.sample_models
        .join(normalizeDOC.out.norm_coverage, by: 0)

    callCNVs(call_cnv_input_ch)

    // Step 10: Collect all CNV calls and filter
    filterCLAMMSCNVs(callCNVs.out.cnvs.map { sample_id, cnv_file -> cnv_file }.collect())

    // Step 11: Convert filtered BED to VCF
    convertClammsToVcf(
        filterCLAMMSCNVs.out.filtered_cnvs,
        sample_file_ch,
        fai_ch
    )

    // Step 12: Sort, compress, and index the VCF
    BGZIP_SORT_INDEX_VCF(convertClammsToVcf.out.vcfs)

    emit:
    sorted_vcf       = BGZIP_SORT_INDEX_VCF.out.sorted_vcf
    sorted_vcf_index = BGZIP_SORT_INDEX_VCF.out.sorted_vcf_index
}
