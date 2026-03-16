#!/usr/bin/env nextflow
nextflow.enable.dsl=2

// =====================================================================================
// MODULE INCLUDES
// =====================================================================================
include { run_Fetch; run_Aggregate; run_Score; run_Database; run_Annotate; run_DenovoTrio; run_DenovoMom; run_DenovoDad; filterINDELIBLE } from './modules/modules-indelible.nf'
include { genReadCounts; calcGC_CANOES; runCANOES; filterCANOESCNVs; convertCanoesToVcf } from './modules/modules-canoes.nf'
include { groupBAMs; gatkDOC; combineDOC; calcGC_XHMM; filterSamples; runPCA; normalisePCA; filterZScore; filterRD; discoverCNVs; genotypeCNVs; splitVCF; filterXHMMCNVs } from './modules/modules-xhmm.nf'
include { generateWindows; samtoolsDOC; normalizeDOC; createPCAData; getPicardQCMetrics; getPicardMeanInsertSize; combinePicardQCMetrics; createCustomRefPanel; trainModels; callCNVs; filterCLAMMSCNVs; convertClammsToVcf } from './modules/modules-clamms.nf'
include { uploadCramFiles; getStaticFiles; checkFileStatus; startAnalysisBatch; checkAnalysisStatus; downloadAnalysisOutput; deleteData; addDragenToolAnnotation } from './modules/modules-icav2-dragen.nf'
include { GENERATE_ACCESS; AUTOBIN; COVERAGE; CREATE_POOLED_REFERENCE; CALL_CNV; EXPORT_RESULTS } from './modules/modules-cnvkit.nf'
include { GENERATE_PLOIDY_PRIORS; PREPROCESS_INTERVALS; ANNOTATE_INTERVALS; COLLECT_READ_COUNTS; FILTER_INTERVALS; DETERMINE_PLOIDY_COHORT; SCATTER_INTERVALS; GERMLINE_CNV_CALLER_COHORT; POSTPROCESS_CALLS } from './modules/modules-gcnv.nf'
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
        calcGC_CANOES(chroms)
        
        // Extract the list of samples from the BAMs directly for the Python script
        bams.map { it -> it[0] + '\n' }.collectFile(name: 'sample_list.txt').set { sample_list }
        
        bams.collectFile() { item -> [ 'bam_list_unsorted.txt', "${item[1]}" + '\n' ] }.set { bam_list }
        genReadCounts(bam_list, chroms)
        genReadCounts.out.chr_reads_cov.join(calcGC_CANOES.out.chr_gc_content).set { chr_canoes_input }
        runCANOES(chr_canoes_input)
        runCANOES.out.chr_cnvs_pass.map { it -> it[1] }.collect().set { all_cnvs_pass }
        filterCANOESCNVs(all_cnvs_pass)
        
        convertCanoesToVcf(filterCANOESCNVs.out.filtered_cnvs, sample_list, fai)
}

workflow RUN_XHMM {
    take: 
        bams
    main:
        groupBAMs(bams.collectFile() { item -> [ 'bam_list_unsorted.txt', "${item[1]}" + '\n' ] })
        gatkDOC(groupBAMs.out.bam_groups.flatMap().map { it -> [it.name[0..-6], it] })
        combineDOC(gatkDOC.out.bam_group_doc.collect { it -> it[1] })
        calcGC_XHMM()
        filterSamples(combineDOC.out.combined_doc, calcGC_XHMM.out.extreme_gc_targets)
        runPCA(filterSamples.out.filtered_centered)
        normalisePCA(filterSamples.out.filtered_centered, runPCA.out.pca_data)
        filterZScore(normalisePCA.out.data_pca_norm)
        filterRD(combineDOC.out.combined_doc, filterSamples.out.excluded_filtered_targets, filterZScore.out.excluded_zscore_targets, filterSamples.out.excluded_filtered_samples, filterZScore.out.excluded_zscore_samples)
        discoverCNVs(filterRD.out.orig_filtered, filterZScore.out.pca_norm_zscore)
        genotypeCNVs(filterRD.out.orig_filtered, filterZScore.out.pca_norm_zscore, discoverCNVs.out.cnvs)
        splitVCF(genotypeCNVs.out.genotypes.flatten().filter { it.name == 'DATA.vcf' })
        filterXHMMCNVs(splitVCF.out.individual_vcfs)
}

workflow RUN_CLAMMS {
    take: 
        bams
        fai
    main:
        generateWindows()
        samtoolsDOC(bams, generateWindows.out.windows)
        normalizeDOC(samtoolsDOC.out.coverage, generateWindows.out.windows)
        normalizeDOC.out.norm_coverage.map { it -> it[1] }.collect().set { norm_coverage_files }
        createPCAData(norm_coverage_files)
        getPicardQCMetrics(bams)
        getPicardMeanInsertSize(bams)
        getPicardQCMetrics.out.qc_metrics.join(getPicardMeanInsertSize.out.ins_size_metrics, by: 0).map { it -> [ it[1], it[2] ] }.flatten().collect().set { picard_metrics }
        combinePicardQCMetrics(picard_metrics)
        createCustomRefPanel(norm_coverage_files, createPCAData.out.pca_data, combinePicardQCMetrics.out.qcs_metrics)
        createCustomRefPanel.out.ref_panel.flatten().map { it -> [ "${it.baseName.replaceAll('.ref.panel.files','')}", it ] }.set { for_training }
        trainModels(for_training, generateWindows.out.windows, norm_coverage_files)
        trainModels.out.sample_models.join(normalizeDOC.out.norm_coverage).set { cllin }
        callCNVs(cllin)
        filterCLAMMSCNVs(callCNVs.out.cnvs.map { it -> it[1] }.collect())
        
        // Extract the list of samples from the BAMs directly to generate the sample list text file
        bams.map { it -> it[0] + '\n' }.collectFile(name: 'sample_list.txt').set { sample_list }
        
        // Pass the output of filterCLAMMSCNVs (assuming it outputs the BED), the sample list, and the FAI index
        convertClammsToVcf(filterCLAMMSCNVs.out.filtered_cnvs, sample_list, fai)
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
        access_ch = GENERATE_ACCESS(fasta)
        first_bam_ch = bams.first().map { it[1] }
        bins_ch = AUTOBIN(fasta, targets, access_ch, refflat, first_bam_ch)
        cov_ch = COVERAGE(bams, bins_ch.target_bed.collect(), bins_ch.antitarget_bed.collect())
        all_covs = cov_ch.target_cov.map { it[1] }.concat(cov_ch.antitarget_cov.map { it[1] }).collect()
        ref_cnn = CREATE_POOLED_REFERENCE(fasta, all_covs)
        cov_pairs = cov_ch.target_cov.join(cov_ch.antitarget_cov)
        calls_ch = CALL_CNV(cov_pairs, ref_cnn.collect())
        EXPORT_RESULTS(calls_ch.results)
}

workflow RUN_GCNV {
    take: 
        bams
        fasta
        fai
        dict
        targets
    main:
        GENERATE_PLOIDY_PRIORS(fai)
        PREPROCESS_INTERVALS(fasta, fai, dict, targets)
        COLLECT_READ_COUNTS(bams, PREPROCESS_INTERVALS.out.interval_list, fasta, fai, dict)
        ANNOTATE_INTERVALS(PREPROCESS_INTERVALS.out.interval_list, fasta, fai, dict)
        
        ch_all_counts = COLLECT_READ_COUNTS.out.counts.map { it[1] }.collect()
        
        FILTER_INTERVALS(PREPROCESS_INTERVALS.out.interval_list, ANNOTATE_INTERVALS.out.annotated_intervals, ch_all_counts)
        DETERMINE_PLOIDY_COHORT(FILTER_INTERVALS.out.filtered_intervals, ch_all_counts, GENERATE_PLOIDY_PRIORS.out.priors)
        SCATTER_INTERVALS(FILTER_INTERVALS.out.filtered_intervals)
        
        ch_scatters = SCATTER_INTERVALS.out.shards.flatten()
        GERMLINE_CNV_CALLER_COHORT(ch_scatters, ANNOTATE_INTERVALS.out.annotated_intervals, ch_all_counts, DETERMINE_PLOIDY_COHORT.out.ploidy_calls)
        
        ch_ploidy_calls = DETERMINE_PLOIDY_COHORT.out.ploidy_calls
        ch_model_shards = GERMLINE_CNV_CALLER_COHORT.out.model.collect()
        ch_call_shards = GERMLINE_CNV_CALLER_COHORT.out.calls.collect()
        
        ch_sample_metadata = COLLECT_READ_COUNTS.out.counts.toSortedList({ a, b -> a[0] <=> b[0] }).flatten().map { it -> it[0] }.indexed()

        POSTPROCESS_CALLS(ch_sample_metadata, ch_model_shards, ch_call_shards, ch_ploidy_calls, dict)
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
