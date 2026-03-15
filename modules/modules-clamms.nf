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
process convertClammsToVCF {
    tag "convert"
    label 'clamms|bedtools'
    publishDir "${outdir}/out_CLAMMS/vcfs", mode: 'copy', overwrite: true

    input:
    path input_file
    path sample_file
    path fai_file
    path log_file

    output:
    path("combined_CLAMMS_output.vcf"), emit: vcf_file

    script:
    """
    python3 clamms_bed_to_vcf.py --input_file ${input_file} --sample_file ${sample_file} --output_dir . --fai_file ${fai_file} --log_file ${log_file} --combined
    echo "Conversion from BED to VCF completed"
    """
}

// Process to bgzip, sort, and index the produced VCF file
process run_BCFtools {
    tag "bcftools"
    label 'clamms|bedtools'
    publishDir "${outdir}/out_CLAMMS/vcfs", mode: 'copy', overwrite: true

    input:
    path vcf_file

    output:
    path("*.sorted.vcf.gz"), emit: sorted_vcf
    path("*.sorted.vcf.gz.tbi"), emit: sorted_vcf_index

    script:
    def sorted_gz = "${vcf_file.simpleName}.sorted.vcf.gz"
    """
    # bgzip the VCF
    bgzip -c ${vcf_file} > ${vcf_file}.gz

    # Sort the bgzipped VCF
    bcftools sort ${vcf_file}.gz -o ${sorted_gz} -O z

    # Index the sorted VCF
    tabix -p vcf ${sorted_gz}

    # Remove the intermediate unsorted gz
    rm ${vcf_file}.gz
    """
}
