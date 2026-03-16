#!/usr/bin/env nextflow
nextflow.enable.dsl=2

// =====================================================================================
// MODULE INCLUDES
// =====================================================================================
include { run_Fetch; run_Aggregate; run_Score; run_Database; run_Annotate; run_DenovoTrio; run_DenovoMom; run_DenovoDad; filterINDELIBLE } from './modules/modules-indelible.nf'
include { CANOES } from './modules/modules-canoes.nf'
include { XHMM } from './modules/modules-xhmm.nf'
include { CLAMMS } from './modules/modules-clamms.nf'
include { uploadCramFiles; getStaticFiles; checkFileStatus; startAnalysisBatch; checkAnalysisStatus; downloadAnalysisOutput; deleteData; addDragenToolAnnotation } from './modules/modules-icav2-dragen.nf'
include { CNVKIT } from './modules/modules-cnvkit.nf'
include { GATK_GCNV } from './modules/modules-gatk-gcnv.nf'
include { runSurvivorMerge } from './modules/modules-survivor.nf'
include { runTruvariCollapse } from './modules/modules-truvari.nf'

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
        run_Fetch(crams)
        run_Aggregate(run_Fetch.out.sc_reads)
        run_Score(run_Aggregate.out.counts)
        run_Database(run_Score.out.database_in.collect())
        run_Annotate(run_Database.out.indel_database, run_Score.out.scores)
        run_DenovoTrio(cram_trios.join(run_Annotate.out.annotated))
        run_DenovoMom(cram_mom.join(run_Annotate.out.annotated))
        run_DenovoDad(cram_dad.join(run_Annotate.out.annotated))
        filterINDELIBLE(run_Annotate.out.annotated)
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
        uploadCramFiles(cramPairs)
        allUploads = uploadCramFiles.out.dataFile.collectFile(name: 'all_data.txt', keepHeader: false, newLine: true)
        getStaticFiles(allUploads)
        checkFileStatus(getStaticFiles.out.dataFile)
        startAnalysisBatch(checkFileStatus.out.dataFile)
        checkAnalysisStatus(startAnalysisBatch.out.dataFile)
        downloadAnalysisOutput(checkAnalysisStatus.out.dataFile)
        deleteData(downloadAnalysisOutput.out.dataFile)
        addDragenToolAnnotation(downloadAnalysisOutput.out.dataFile)
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
            
        runSurvivorMerge(grouped_vcfs)
}

workflow RUN_TRUVARI {
    take: 
        vcf_ch
    main:
        // Group by sample_id, then enforce that a sample MUST have VCFs from >= 2 different callers to proceed
        grouped_vcfs = vcf_ch
            .groupTuple()
            .filter { sample_id, vcfs -> vcfs.size() >= 2 }
            
        runTruvariCollapse(grouped_vcfs)
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
            Channel.fromFilePairs("${params.cramFilePairsUploadPath}", checkIfExists: true) { file -> file.name.replaceAll(/.cram|.crai|.cram.crai$/, '') }
                .set { ch_cramPairs }
            RUN_DRAGEN(ch_cramPairs)
            break

        case['cnvkit']:
            Channel.fromPath(params.bams.replace(".bam", ".{bam,bai}"))
                .map { it -> [it.baseName, it] }.groupTuple(size: 2)
                .map { id, files -> [id, files.find { it.extension == 'bam' }, files.find { it.extension == 'bai' }] }
                .filter { id, bam, bai -> 
                    (params.test_list != "all" && params.test_list != "") ? params.test_list.tokenize(',').contains(id) : true 
                }
                .set { ch_bams_filtered }
            
            ch_bams_final = (params.test_size.toInteger() > 0) ? ch_bams_filtered.take(params.test_size.toInteger()) : ch_bams_filtered
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
"""
            break
    }
}
