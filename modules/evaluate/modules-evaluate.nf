#!/usr/bin/env nextflow
nextflow.enable.dsl=2

// ================================================================================
//  EVALUATE MODULE
//  Assess CNV caller performance (precision and sensitivity) at the probe level.
//
//  Pipeline:
//    1. VCF_TO_BED   – convert each per-sample VCF to a 5-column BED file
//    2. COMBINE_BEDS – concatenate all per-sample BED files into a single
//                      unified call-set BED
//    3. EVALUATE_CALLER – compare the unified call set against the truth set
//                         using evaluate_caller_performance.py and report
//                         precision / sensitivity metrics
//
//  BED column order expected by steps 2 & 3 (0-based, half-open):
//    CHR  START  STOP  CNV_TYPE  SAMPLE_ID
// ================================================================================

outdir = file(params.outdir, type: 'dir')

// ---------------------------------------------------------------------------------
// Process 1: Convert a single per-sample VCF to a 5-column BED file
// ---------------------------------------------------------------------------------
process VCF_TO_BED {
    tag "${vcf.simpleName}"
    label 'pysam'
    publishDir "${outdir}/out_EVALUATE", mode: 'copy', overwrite: true

    input:
    path vcf

    output:
    path("*.bed"), emit: bed

    script:
    // Strip common suffixes to obtain a clean sample name
    def sample_name = vcf.name
        .replaceAll(/\.(sorted|annotated|normalised)(\.vcf(\.gz)?)?$/, '')
        .replaceAll(/\.vcf(\.gz)?$/, '')
    """
    vcf_to_bed.py --vcf ${vcf} --output ${sample_name}.bed
    """
}

// ---------------------------------------------------------------------------------
// Process 2: Combine all per-sample BED files into one unified call-set BED
// ---------------------------------------------------------------------------------
process COMBINE_BEDS {
    tag "${caller_name}"
    label 'pysam'
    publishDir "${outdir}/out_EVALUATE", mode: 'copy', overwrite: true

    input:
    path beds
    val  caller_name

    output:
    path("${caller_name}_callset.bed"), emit: combined_bed

    script:
    """
    cat ${beds} > ${caller_name}_callset.bed
    """
}

// ---------------------------------------------------------------------------------
// Process 3: Run evaluate_caller_performance.py and write metrics to a text file
// ---------------------------------------------------------------------------------
process EVALUATE_CALLER {
    tag "${caller_name}"
    label 'pysam'
    publishDir "${outdir}/out_EVALUATE", mode: 'copy', overwrite: true

    input:
    path callset_bed
    path truth_bed
    path probes_bed
    val  caller_name

    output:
    path("${caller_name}_performance.txt"), emit: report

    script:
    """
    evaluate_caller_performance.py \\
        --truth_bed   ${truth_bed}   \\
        --callset_bed ${callset_bed} \\
        --probes_bed  ${probes_bed}  \\
        --output      ${caller_name}_performance.txt
    """
}

// ---------------------------------------------------------------------------------
// EVALUATE workflow
// ---------------------------------------------------------------------------------
workflow EVALUATE {
    take:
    vcf_ch      // channel of per-sample VCF paths (sorted, annotated, or normalised)
    truth_bed   // path: truth set BED (CHR, START, STOP, CNV_TYPE, SAMPLE_ID)
    probes_bed  // path: capture target BED (CHR, START, STOP[, ...])
    caller_name // val:  caller label used in output filenames (e.g. CANOES)

    main:
    vcf_ch.ifEmpty { error "EVALUATE received no VCF files. Check that vcf_ch contains VCF file paths." }
    VCF_TO_BED(vcf_ch)
    COMBINE_BEDS(VCF_TO_BED.out.bed.collect(), caller_name)
    EVALUATE_CALLER(
        COMBINE_BEDS.out.combined_bed,
        truth_bed,
        probes_bed,
        caller_name
    )

    emit:
    performance_report = EVALUATE_CALLER.out.report
}
