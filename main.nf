#!/usr/bin/env nextflow
nextflow.enable.dsl=2

// =====================================================================================
// MODULE INCLUDES
// =====================================================================================
include { INDELIBLE } from './modules/modules-indelible.nf'
include { CANOES } from './modules/modules-canoes.nf'
include { XHMM } from './modules/modules-xhmm.nf'
include { CLAMMS } from './modules/modules-clamms.nf'
include { DRAGEN } from './modules/modules-icav2-dragen.nf'
include { CNVKIT } from './modules/modules-cnvkit.nf'
include { GATK_GCNV } from './modules/modules-gatk-gcnv.nf'
include { SURVIVOR } from './modules/modules-survivor.nf'
include { TRUVARI } from './modules/modules-truvari.nf'
include { FEATURE_EXTRACTION } from './modules/modules-feature-extraction.nf'
include { NORMALISE } from './modules/modules-normalise.nf'
include { TRAIN } from './modules/modules-train.nf'
include { EVALUATE } from './modules/modules-evaluate.nf'

// =====================================================================================
// GLOBAL SETUP
// =====================================================================================
params.outdir = './output' // Change this to your desired output path
outdir = file(params.outdir, type: 'dir')
outdir.mkdir()
workflow_mode = params.workflow

// =====================================================================================
// HELPER FUNCTION: GATHER VCFS FOR CONSENSUS MODULES
// =====================================================================================
def gather_vcfs() {
    def ch = Channel.empty()
    def dir_count = 0
    
    // Helper closure to reliably extract the sample ID from various caller filenames.
    // Also handles DRAGEN Germline Enrichment outputs named ${sample_id}.cnv.vcf.gz.
    def get_id = { f -> f.name.replaceAll(/_(CANOES|CLAMMS|XHMM|CNVKIT|GCNV|DRAGEN|INDELIBLE).*/, '').replaceAll(/\.cnv\.vcf(\.gz)?$/i, '').replaceAll(/\.vcf(\.gz)?$/i, '') }

    if (params.get('canoes_dir', false))    { ch = ch.mix(Channel.fromPath(params.canoes_dir + "/*.vcf*").map { f -> [get_id(f), f] }); dir_count++ }
    if (params.get('clamms_dir', false))    { ch = ch.mix(Channel.fromPath(params.clamms_dir + "/*.vcf*").map { f -> [get_id(f), f] }); dir_count++ }
    if (params.get('xhmm_dir', false))      { ch = ch.mix(Channel.fromPath(params.xhmm_dir + "/*.vcf*").map   { f -> [get_id(f), f] }); dir_count++ }
    if (params.get('cnvkit_dir', false))    { ch = ch.mix(Channel.fromPath(params.cnvkit_dir + "/*.vcf*").map { f -> [get_id(f), f] }); dir_count++ }
    if (params.get('gcnv_dir', false))      { ch = ch.mix(Channel.fromPath(params.gcnv_dir + "/*.vcf*").map   { f -> [get_id(f), f] }); dir_count++ }
    if (params.get('dragen_dir', false))    { ch = ch.mix(Channel.fromPath(params.dragen_dir + "/*.vcf*").map { f -> [get_id(f), f] }); dir_count++ }
    if (params.get('indelible_dir', false)) { ch = ch.mix(Channel.fromPath(params.indelible_dir + "/*.vcf*").map { f -> [get_id(f), f] }); dir_count++ }
    
    if (dir_count < 2) {
        exit 1, "Error: You must provide VCF directories for at least TWO different callers (e.g., --canoes_dir and --clamms_dir) to run consensus modules."
    }
    
    return ch
}

// =====================================================================================
// SUB-WORKFLOW DEFINITIONS
// =====================================================================================

workflow RUN_INDELIBLE {
    take: 
        crams
        cram_trios
        cram_mom
        cram_dad
    main:
        INDELIBLE(crams, cram_trios, cram_mom, cram_dad)
        if (params.get('truth_bed', false) && params.get('probes_bed', false)) {
            EVALUATE(
                INDELIBLE.out.sorted_vcf.flatten(),
                file(params.truth_bed),
                file(params.probes_bed),
                'INDELIBLE'
            )
        }
}

workflow RUN_CANOES {
    take: 
        bams
        chroms
        fai
    main:
        bam_list = bams.collectFile() { item -> [ 'bam_list_unsorted.txt', "${item[1]}" + '\n' ] }
        CANOES(bam_list, fai, Channel.from(chroms))
        if (params.get('truth_bed', false) && params.get('probes_bed', false)) {
            EVALUATE(
                CANOES.out.sorted_vcf.flatten(),
                file(params.truth_bed),
                file(params.probes_bed),
                'CANOES'
            )
        }
}

workflow RUN_XHMM {
    take: 
        bams
    main:
        bam_list = bams.collectFile() { item -> [ 'bam_list_unsorted.txt', "${item[1]}" + '\n' ] }
        XHMM(bam_list)
        if (params.get('truth_bed', false) && params.get('probes_bed', false)) {
            EVALUATE(
                XHMM.out.sorted_vcf.flatten(),
                file(params.truth_bed),
                file(params.probes_bed),
                'XHMM'
            )
        }
}

workflow RUN_CLAMMS {
    take: 
        bams
        fai
    main:
        sample_list = bams.map { it -> it[0] + '\n' }.collectFile(name: 'sample_list.txt')
        CLAMMS(bams, fai, sample_list)
        if (params.get('truth_bed', false) && params.get('probes_bed', false)) {
            EVALUATE(
                CLAMMS.out.sorted_vcf.flatten(),
                file(params.truth_bed),
                file(params.probes_bed),
                'CLAMMS'
            )
        }
}

workflow RUN_DRAGEN {
    take: 
        cramPairs
    main:
        DRAGEN(cramPairs)
        if (params.get('truth_bed', false) && params.get('probes_bed', false)) {
            EVALUATE(
                DRAGEN.out.sorted_vcf.flatten(),
                file(params.truth_bed),
                file(params.probes_bed),
                'DRAGEN'
            )
        }
}

workflow RUN_CNVKIT {
    take: 
        bams
        fasta
        targets
        refflat
    main:
        CNVKIT(bams, fasta, targets, refflat, bams.first().map { it[1] })
        if (params.get('truth_bed', false) && params.get('probes_bed', false)) {
            EVALUATE(
                CNVKIT.out.sorted_vcf.flatten(),
                file(params.truth_bed),
                file(params.probes_bed),
                'CNVKIT'
            )
        }
}

workflow RUN_GCNV {
    take: 
        bams
        fasta
        fai
        dict
        targets
    main:
        GATK_GCNV(bams, fasta, fai, dict, targets)
        if (params.get('truth_bed', false) && params.get('probes_bed', false)) {
            EVALUATE(
                GATK_GCNV.out.sorted_vcf.flatten(),
                file(params.truth_bed),
                file(params.probes_bed),
                'GCNV'
            )
        }
}

workflow RUN_EVALUATE {
    take:
        vcf_ch
        truth_bed
        probes_bed
        caller_name
    main:
        EVALUATE(vcf_ch, truth_bed, probes_bed, caller_name)
}

workflow RUN_SURVIVOR {
    take: 
        vcf_ch
    main:
        // Group by sample_id, then enforce that a sample MUST have VCFs from >= 2 different callers to proceed
        grouped_vcfs = vcf_ch
            .groupTuple()
            .filter { sample_id, vcfs -> vcfs.size() >= 2 }
            
        SURVIVOR(grouped_vcfs)
}

workflow RUN_TRUVARI {
    take: 
        vcf_ch
    main:
        // Group by sample_id, then enforce that a sample MUST have VCFs from >= 2 different callers to proceed
        grouped_vcfs = vcf_ch
            .groupTuple()
            .filter { sample_id, vcfs -> vcfs.size() >= 2 }
            
        TRUVARI(grouped_vcfs)
}

workflow RUN_FEATURE_EXTRACTION {
    take:
        // Channel of tuples:
        // [ sample_id, merged_vcf, collapsed_vcf, tool_vcfs_str, merger_mode,
        //   bam_file, reference_fasta, bed_file, mappability_file, indelible_counts ]
        // collapsed_vcf should be [] when not available (SURVIVOR mode or absent).
        feature_inputs_ch
    main:
        FEATURE_EXTRACTION(feature_inputs_ch)
}

// Helper closure: build the feature_inputs tuple for a single merged VCF.
// collapsed_vcf_f: a file path or [] (empty list)
def build_feature_inputs(sample_id, merged_vcf, collapsed_vcf_f) {
    def parts = []
    if (params.get('canoes_norm_dir',    false)) parts << "canoes=${params.canoes_norm_dir}/${sample_id}_CANOES.normalised.vcf.gz"
    if (params.get('clamms_norm_dir',    false)) parts << "clamms=${params.clamms_norm_dir}/${sample_id}_CLAMMS.normalised.vcf.gz"
    if (params.get('xhmm_norm_dir',      false)) parts << "xhmm=${params.xhmm_norm_dir}/${sample_id}_XHMM.normalised.vcf.gz"
    if (params.get('cnvkit_norm_dir',    false)) parts << "cnvkit=${params.cnvkit_norm_dir}/${sample_id}_CNVKIT.normalised.vcf.gz"
    if (params.get('gcnv_norm_dir',      false)) parts << "gatk_gcnv=${params.gcnv_norm_dir}/${sample_id}_GCNV.normalised.vcf.gz"
    if (params.get('dragen_norm_dir',    false)) parts << "dragen=${params.dragen_norm_dir}/${sample_id}_DRAGEN.normalised.vcf.gz"
    if (params.get('indelible_norm_dir', false)) parts << "indelible=${params.indelible_norm_dir}/${sample_id}_INDELIBLE.normalised.vcf.gz"
    def tool_vcfs_str = parts.join(',')

    def bam_f       = params.get('bam_file',         false) ? file(params.bam_file)         : []
    def fasta_f     = params.get('reference_fasta',  false) ? file(params.reference_fasta)  : []
    def bed_f       = params.get('bed_file',         false) ? file(params.bed_file)         : []
    def map_f       = params.get('mappability_file', false) ? file(params.mappability_file) : []
    def indelible_f = params.get('indelible_counts', false) ? file(params.indelible_counts) : []
    def mode        = params.get('merger_mode', 'survivor')

    return [sample_id, merged_vcf, collapsed_vcf_f, tool_vcfs_str, mode,
            bam_f, fasta_f, bed_f, map_f, indelible_f]
}

workflow RUN_NORMALISE {
    take:
        // Channel of tuples: [ sample_name, vcf_path, caller_name ]
        vcf_ch
    main:
        NORMALISE(vcf_ch)
}

workflow RUN_TRAIN {
    take:
        // Channel of individual *_features.tsv file paths (one per sample).
        features_tsv_ch
        // Single-element channel containing the path to the truth-labels TSV.
        // Required columns: sample_id, chrom, start, end, truth_label
        truth_labels_ch
    main:
        TRAIN(features_tsv_ch, truth_labels_ch)
}

workflow RUN_SURVIVOR_WITH_FEATURES {
    take:
        vcf_ch
    main:
        grouped_vcfs = vcf_ch
            .groupTuple()
            .filter { sample_id, vcfs -> vcfs.size() >= 2 }
        SURVIVOR(grouped_vcfs)
        feature_inputs_ch = SURVIVOR.out.union_vcf
            .map { sample_id, merged_vcf -> build_feature_inputs(sample_id, merged_vcf, []) }
        FEATURE_EXTRACTION(feature_inputs_ch)
}

workflow RUN_TRUVARI_WITH_FEATURES {
    take:
        vcf_ch
    main:
        grouped_vcfs = vcf_ch
            .groupTuple()
            .filter { sample_id, vcfs -> vcfs.size() >= 2 }
        TRUVARI(grouped_vcfs)
        feature_inputs_ch = TRUVARI.out.merged_vcf
            .join(TRUVARI.out.collapsed_vcf)
            .map { sample_id, merged_vcf, collapsed_vcf -> build_feature_inputs(sample_id, merged_vcf, collapsed_vcf) }
        FEATURE_EXTRACTION(feature_inputs_ch)
}

// =====================================================================================
// WORKFLOW SWITCH - FULLY INDEPENDENT EXECUTION BLOCKS
// =====================================================================================

workflow {
    switch (workflow_mode) {
        
        case['indelible']:
            Channel.fromFilePairs([params.crams + '/*{.cram,.cram.crai}'])
                .map { it -> [ it[0][0..-6], it[1][0], it[1][1] ] }
                .filter { it -> it[1] =~ '_01_1' }
                .set { ch_crams }

            Channel.fromFilePairs([params.crams + '/*{_01_1,_02_2,_03_3}*'], size: 6)
                .map { it -> [ it[0], it[1][0], it[1][1], it[1][2], it[1][3], it[1][4], it[1][5] ] }
                .set { ch_cram_trios }

            Channel.fromFilePairs([params.crams + '/*{_01_1,_02_2,_03_3}*'], size: -1)
                .filter { it -> it[1].size() == 4 && it[1][3] =~ '02_2.cram' }
                .map { it -> [ it[0], it[1][0], it[1][1], it[1][2], it[1][3] ] }
                .set { ch_cram_mom }

            Channel.fromFilePairs([params.crams + '/*{_01_1,_02_2,_03_3}*'], size: -1)
                .filter { it -> it[1].size() == 4 && it[1][3] =~ '03_3.cram' }
                .map { it -> [ it[0], it[1][0], it[1][1], it[1][2], it[1][3] ] }
                .set { ch_cram_dad }

            RUN_INDELIBLE(ch_crams, ch_cram_trios, ch_cram_mom, ch_cram_dad)
            break

        case['canoes']:
            chroms = (1..22).toList().collect { 'chr' + "${it}" }
            Channel.fromPath(params.samplesheet_bams)
                .splitCsv(header: true, sep: '\t')
                .map { row -> [ "${row.SampleID}", "${row.BAM}", "${row.BAM}".replaceAll("\\b.bam\\b",".bam.bai") ] }
                .set { ch_bams }
            RUN_CANOES(ch_bams, chroms, file(params.fai))
            break

        case['xhmm']:
            Channel.fromPath(params.samplesheet_bams)
                .splitCsv(header: true, sep: '\t')
                .map { row -> [ "${row.SampleID}", "${row.BAM}", "${row.BAM}".replaceAll("\\b.bam\\b",".bam.bai") ] }
                .set { ch_bams }
            RUN_XHMM(ch_bams)
            break

        case['clamms']:
            Channel.fromPath(params.samplesheet_bams)
                .splitCsv(header: true, sep: '\t')
                .map { row -> [ "${row.SampleID}", "${row.BAM}", "${row.BAM}".replaceAll("\\b.bam\\b",".bam.bai") ] }
                .set { ch_bams }
            RUN_CLAMMS(ch_bams, file(params.fai))
            break

        case['dragen']:
            // Support both CRAM/CRAI and BAM/BAI file pairs
            Channel.fromFilePairs(
                    ["${params.cramFilePairsUploadPath}"],
                    checkIfExists: true
                ) { file -> file.name.replaceAll(/\.(cram\.crai|bam\.bai|cram|crai|bam|bai)$/, '') }
                .set { ch_cramPairs }
            RUN_DRAGEN(ch_cramPairs)
            break

        case['cnvkit']:
            Channel.fromPath(params.bams.replace(".bam", ".{bam,bai}"))
                .map { it -> [it.baseName, it] }.groupTuple(size: 2)
                .map { id, files -> [id, files.find { it.extension == 'bam' }, files.find { it.extension == 'bai' }] }
                .filter { id, bam, bai -> 
                    def tl = params.get('test_list', 'all')
                    (tl != "all" && tl != "") ? tl.tokenize(',').contains(id) : true 
                }
                .set { ch_bams_filtered }
            
            ch_bams_final = (params.get('test_size', -1).toInteger() > 0) ? ch_bams_filtered.take(params.get('test_size', -1).toInteger()) : ch_bams_filtered
            RUN_CNVKIT(ch_bams_final, file(params.fasta), file(params.targets), file(params.refflat))
            break

        case['gcnv']:
            Channel.fromPath(params.samples_path)
                .map { it -> 
                    def index = it.name.endsWith('.bam') ? "${it}.bai" : "${it}.crai"
                    return [ it.baseName, it, file(index) ] 
                }
                .set { ch_bams }
            RUN_GCNV(ch_bams, file(params.fasta), file(params.fai), file(params.dict), file(params.exome_targets))
            break

        case['survivor']:
            def ch_vcfs = gather_vcfs()
            RUN_SURVIVOR(ch_vcfs)
            break

        case['truvari']:
            def ch_vcfs = gather_vcfs()
            RUN_TRUVARI(ch_vcfs)
            break

        case['survivor_with_features']:
            def ch_vcfs = gather_vcfs()
            RUN_SURVIVOR_WITH_FEATURES(ch_vcfs)
            break

        case['truvari_with_features']:
            def ch_vcfs = gather_vcfs()
            RUN_TRUVARI_WITH_FEATURES(ch_vcfs)
            break

        case['normalise']:
            // Normalise quality scores for pre-existing caller VCFs without running
            // the full caller pipeline.
            // Required: --vcf_dir  (directory containing raw .vcf / .vcf.gz files)
            //           --caller   (CANOES | CLAMMS | XHMM | GATK | CNVKIT | DRAGEN | INDELIBLE)
            Channel.fromPath(params.vcf_dir + '/*.vcf*')
                .map { f ->
                    def sample_name = f.name
                        .replaceAll(/_(CANOES|CLAMMS|XHMM|CNVKIT|GCNV|DRAGEN|INDELIBLE).*/, '')
                        .replaceAll(/\.(sorted\.)?vcf(\.gz)?$/, '')
                    [sample_name, f, params.caller.toUpperCase()]
                }
                .set { ch_vcfs }
            RUN_NORMALISE(ch_vcfs)
            break

        case['feature_extraction']:
            // Run feature extraction on an already-produced merged VCF.
            // Required: --merged_vcf_dir  (directory containing <sample_id>_survivor_union.vcf
            //                              or <sample_id>_truvari_merged.vcf files)
            // Optional annotation params (any subset may be omitted):
            //   --canoes_norm_dir, --clamms_norm_dir, --xhmm_norm_dir, --cnvkit_norm_dir,
            //   --gcnv_norm_dir, --dragen_norm_dir, --indelible_norm_dir
            //   --bam_file, --reference_fasta, --bed_file, --mappability_file,
            //   --indelible_counts, --merger_mode (default: survivor)
            if (params.get('merger_mode', 'survivor') == 'truvari') {
                Channel.fromPath(params.merged_vcf_dir + '/**/*_truvari_merged.vcf*')
                    .map { f ->
                        def sample_id = f.name.replaceAll(/_truvari.*/, '').replaceAll(/\.vcf(\.gz)?$/i, '')
                        def collapsed_f = file("${f.parent}/${sample_id}_truvari_collapsed.vcf")
                        def collapsed_vcf_f = collapsed_f.exists() ? collapsed_f : []
                        build_feature_inputs(sample_id, f, collapsed_vcf_f)
                    }
                    .set { ch_feature_inputs }
            } else {
                Channel.fromPath(params.merged_vcf_dir + '/**/*_survivor_union.vcf*')
                    .map { f ->
                        def sample_id = f.name.replaceAll(/_survivor.*/, '').replaceAll(/\.vcf(\.gz)?$/i, '')
                        build_feature_inputs(sample_id, f, [])
                    }
                    .set { ch_feature_inputs }
            }
            RUN_FEATURE_EXTRACTION(ch_feature_inputs)
            break

        case['train']:
            // Train an XGBoost classifier on CNV feature matrices produced by
            // the feature_extraction workflow.
            // Required: --features_dir   (directory containing *_features.tsv files
            //                             produced by the feature_extraction workflow)
            //           --truth_labels   (TSV file with columns:
            //                             sample_id, chrom, start, end, truth_label
            //                             where truth_label=1 means true CNV)
            Channel.fromPath(params.features_dir + '/**/*_features.tsv')
                .set { ch_features }
            Channel.value(file(params.truth_labels))
                .set { ch_truth }
            RUN_TRAIN(ch_features, ch_truth)
            break

        case['evaluate']:
            // Evaluate CNV caller performance (precision and sensitivity) using
            // pre-existing per-sample VCFs produced by one of the 7 callers.
            // VCFs are converted to BED, merged into a unified call set, and
            // compared against the truth set at the probe (capture target) level.
            //
            // Required: --vcf_dir    (directory containing per-sample .vcf / .vcf.gz files)
            //           --caller     (caller label used in output names, e.g. CANOES)
            //           --truth_bed  (truth set BED: CHR, START, STOP, CNV_TYPE, SAMPLE_ID)
            //           --probes_bed (capture target BED: CHR, START, STOP[, ...])
            Channel.fromPath(params.vcf_dir + '/*.vcf*')
                .set { ch_vcfs }
            RUN_EVALUATE(
                ch_vcfs,
                file(params.truth_bed),
                file(params.probes_bed),
                params.caller
            )
            break

        default:
            exit 1, """
OOOPS!! SEEMS LIKE WE HAVE A WORKFLOW ERROR!

No workflow 'mode' given!
Please use one of the following options for workflows:
    --workflow indelible    
    --workflow canoes       
    --workflow xhmm       
    --workflow clamms      
    --workflow dragen      
    --workflow cnvkit      
    --workflow gcnv       
    --workflow survivor     
    --workflow truvari
    --workflow survivor_with_features
    --workflow truvari_with_features
    --workflow feature_extraction
    --workflow normalise
    --workflow train
    --workflow evaluate
"""
            break
    }
}
