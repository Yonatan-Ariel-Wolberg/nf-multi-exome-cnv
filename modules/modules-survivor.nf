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

process runSurvivorMerge {
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
    SURVIVOR merge \\
        \$list_file \\
        1000 \\      # max distance
        1 \\         # min support
        1 \\         # use type
        0 \\         # use strand
        0 \\         # estimate distance
        30 \\        # min size
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
        runSurvivorMerge(grouped_vcfs)

    emit:
        merged_vcf = runSurvivorMerge.out.merged_vcf
}
