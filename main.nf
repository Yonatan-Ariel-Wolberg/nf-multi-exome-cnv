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
    
    // Helper closure to reliably extract the sample ID from various caller filenames
    def get_id = { f -> f.name.replaceAll(/_(CANOES|CLAMMS|XHMM|CNVKIT|GCNV|DRAGEN|INDELIBLE).*/, '').replaceAll(/\.vcf(\.gz)?$/i, '') }

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
}

workflow RUN_CANOES {
    take: 
        bams
        chroms
        fai
    main:
        bam_list = bams.collectFile() { item -> [ 'bam_list_unsorted.txt', "${item[1]}" + '\n' ] }
        CANOES(bam_list, fai, Channel.from(chroms))
}

workflow RUN_XHMM {
    take: 
        bams
    main:
        bam_list = bams.collectFile() { item -> [ 'bam_list_unsorted.txt', "${item[1]}" + '\n' ] }
        XHMM(bam_list)
}

workflow RUN_CLAMMS {
    take: 
        bams
        fai
    main:
        sample_list = bams.map { it -> it[0] + '\n' }.collectFile(name: 'sample_list.txt')
        CLAMMS(bams, fai, sample_list)
}

workflow RUN_DRAGEN {
    take: 
        cramPairs
    main:
        DRAGEN(cramPairs)
}

workflow RUN_CNVKIT {
    take: 
        bams
        fasta
        targets
        refflat
    main:
        CNVKIT(bams, fasta, targets, refflat, bams.first().map { it[1] })
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
        // [ sample_id, merged_vcf ]
        // Optional annotation params are taken from the global params block.
        merged_vcf_ch
    main:
        // Build the tool_vcfs string from any per-caller normalised VCF dirs
        // provided via params.  Only callers whose params are set are included.
        feature_inputs = merged_vcf_ch.map { sample_id, merged_vcf ->
            def parts = []
            if (params.get('canoes_norm_dir',    false)) parts << "canoes=${params.canoes_norm_dir}/${sample_id}_CANOES.normalised.vcf.gz"
            if (params.get('clamms_norm_dir',    false)) parts << "clamms=${params.clamms_norm_dir}/${sample_id}_CLAMMS.normalised.vcf.gz"
            if (params.get('xhmm_norm_dir',      false)) parts << "xhmm=${params.xhmm_norm_dir}/${sample_id}_XHMM.normalised.vcf.gz"
            if (params.get('cnvkit_norm_dir',    false)) parts << "cnvkit=${params.cnvkit_norm_dir}/${sample_id}_CNVKIT.normalised.vcf.gz"
            if (params.get('gcnv_norm_dir',      false)) parts << "gatk_gcnv=${params.gcnv_norm_dir}/${sample_id}_GCNV.normalised.vcf.gz"
            if (params.get('dragen_norm_dir',    false)) parts << "dragen=${params.dragen_norm_dir}/${sample_id}_DRAGEN.normalised.vcf.gz"
            if (params.get('indelible_norm_dir', false)) parts << "indelible=${params.indelible_norm_dir}/${sample_id}_INDELIBLE.normalised.vcf.gz"
            def tool_vcfs_str = parts.join(',')

            def bam_f          = params.get('bam_file',          false) ? file(params.bam_file)          : []
            def fasta_f        = params.get('reference_fasta',   false) ? file(params.reference_fasta)   : []
            def bed_f          = params.get('bed_file',          false) ? file(params.bed_file)          : []
            def map_f          = params.get('mappability_file',  false) ? file(params.mappability_file)  : []
            def indelible_f    = params.get('indelible_counts',  false) ? file(params.indelible_counts)  : []
            def mode           = params.get('merger_mode', 'survivor')

            return [sample_id, merged_vcf, tool_vcfs_str, mode,
                    bam_f, fasta_f, bed_f, map_f, indelible_f]
        }

        FEATURE_EXTRACTION(feature_inputs)
}

workflow RUN_NORMALISE {
    take:
        // Channel of tuples: [ sample_name, vcf_path, caller_name ]
        vcf_ch
    main:
        NORMALISE(vcf_ch)
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
            def vcf_pattern = params.get('merger_mode', 'survivor') == 'truvari'
                ? params.merged_vcf_dir + '/**/*truvari*.vcf*'
                : params.merged_vcf_dir + '/**/*survivor*.vcf*'
            Channel.fromPath(vcf_pattern)
                .map { f -> [f.name.replaceAll(/_survivor.*|_truvari.*/, '').replaceAll(/\.vcf(\.gz)?$/i, ''), f] }
                .set { ch_merged }
            RUN_FEATURE_EXTRACTION(ch_merged)
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
    --workflow feature_extraction
    --workflow normalise
"""
            break
    }
}
