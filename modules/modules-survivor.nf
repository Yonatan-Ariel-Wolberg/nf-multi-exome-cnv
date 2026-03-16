#!/usr/bin/env nextflow
nextflow.enable.dsl=2

// =====================================================================================
// GLOBAL FILE INSTANTIATION
// =====================================================================================
params.outdir = './output' // Change this to your desired output path
outdir = file(params.outdir, type: 'dir')

// =====================================================================================
// PROCESSES FOR SURVIVOR
// =====================================================================================

process RUN_SURVIVOR_MERGE {
    tag "${sample_id}"
    label 'survivor'
    publishDir "${outdir}/out_SURVIVOR/${sample_id}", mode: 'copy', overwrite: true

    input:
    // Takes a sample ID and a list of VCF files (from multiple callers) associated with that sample
    tuple val(sample_id), path(vcfs)

    output:
    tuple val(sample_id), path("${sample_id}_survivor_merged.vcf"), emit: merged_vcf

    script:
    """
    # Create a list file for SURVIVOR
    list_file="${sample_id}_vcfs.list"
    rm -f "\$list_file"  # Remove any existing list file to avoid conflicts

    # Loop through the VCF files to process
    for vcf in ${vcfs}; do
        if [[ "\$vcf" == *.gz ]]; then
            uncompressed="\${vcf%.gz}"
            gunzip -c "\$vcf" > "\$uncompressed"  # Decompress on the fly
            echo "\$uncompressed" >> "\$list_file"
        else
            echo "\$vcf" >> "\$list_file"
        fi
    done

    # Run SURVIVOR merge with the specified parameters
    # Args: list max_dist min_support use_type use_strand est_dist min_sv_size output
    SURVIVOR merge \\
        \$list_file \\
        1000 \\
        1 \\
        1 \\
        0 \\
        0 \\
        30 \\
        "${sample_id}_survivor_merged.vcf"
    """
}

// =====================================================================================
// SUB-WORKFLOW TO CHAIN THE PROCESSES TOGETHER
// =====================================================================================

workflow SURVIVOR {
    take:
        grouped_vcfs

    main:
        RUN_SURVIVOR_MERGE(grouped_vcfs)

    emit:
        merged_vcf = RUN_SURVIVOR_MERGE.out.merged_vcf
}
