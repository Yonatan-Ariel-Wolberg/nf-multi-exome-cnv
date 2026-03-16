#!/usr/bin/env nextflow
nextflow.enable.dsl=2

// =====================================================================================
// GLOBAL FILE INSTANTIATION
// =====================================================================================
params {
    outdir = './output' // Change to your desired output path
    bin_length = 1000 // Example parameter for bin length
    padding = 100 // Example parameter for padding
    scatter_count = 10 // Example parameter for scatter count
    is_wgs = false // Set to false indicating the analysis is exome sequencing
}

// Define output directory
outdir = file(params.outdir, type: 'dir')

// Process to generate ploidy priors
process GENERATE_PLOIDY_PRIORS {
    tag "Ploidy_Priors"
    publishDir "${outdir}/out_GCNV/priors", mode: 'copy', overwrite: true

    input: 
    path fai
    
    output: 
    path "ploidy_priors.tsv", emit: priors
    
    script:
    """
    # Create header with 4 ploidy states (0, 1, 2, 3)
    echo -e "CONTIG_NAME\tPLOIDY_PRIOR_0\tPLOIDY_PRIOR_1\tPLOIDY_PRIOR_2\tPLOIDY_PRIOR_3" > ploidy_priors.tsv

    # Default for all contigs: high priority for Ploidy 2 (Autosomes)
    awk -v OFS="\t" '{print \$1, "0.01", "0.01", "0.97", "0.01"}' $fai >> ploidy_priors.tsv

    # Adjust Sex Chromosomes (X and Y) for mixed cohort
    sed -i 's/^chrX\t0.01\t0.01\t0.97\t0.01/chrX\t0.01\t0.49\t0.49\t0.01/' ploidy_priors.tsv
    sed -i 's/^chrY\t0.01\t0.01\t0.97\t0.01/chrY\t0.50\t0.50\t0.00\t0.00/' ploidy_priors.tsv
    """
}

// Process to preprocess intervals
process PREPROCESS_INTERVALS {
    tag "Preprocess_Intervals"
    label 'gatk'
    
    input: 
    path fasta
    path fai
    path dict
    path targets
    
    output: 
    path "preprocessed.interval_list", emit: interval_list
    
    script:
    """
    gatk --java-options "-Xmx4g" PreprocessIntervals \\
        -R $fasta \\
        -L $targets \\
        --bin-length ${params.bin_length} \\
        --padding ${params.padding} \\
        -imr OVERLAPPING_ONLY \\
        -O preprocessed.interval_list
    """
}

// Process to annotate intervals
process ANNOTATE_INTERVALS {
    tag "Annotate_Intervals"
    label 'gatk'
    
    input: 
    path intervals
    path fasta
    path fai
    path dict
    
    output: 
    path "annotated.tsv", emit: annotated_intervals
    
    script:
    """
    gatk --java-options "-Xmx4g" AnnotateIntervals \\
        -L $intervals -R $fasta -imr OVERLAPPING_ONLY -O annotated.tsv
    """
}

// Process to collect read counts
process COLLECT_READ_COUNTS {
    tag "$sample_id"
    label 'gatk'
    
    input: 
    tuple val(sample_id), path(reads), path(index)
    path intervals
    path fasta
    path fai
    path dict

    output: 
    tuple val(sample_id), path("${sample_id}.counts.hdf5"), emit: counts

    script:
    """
    gatk --java-options "-Xmx4g" CollectReadCounts \\
        -L $intervals -R $fasta -I $reads \\
        --format HDF5 -imr OVERLAPPING_ONLY \\
        -O ${sample_id}.counts.hdf5
    """
}

// Process to filter intervals
process FILTER_INTERVALS {
    tag "Filter_Intervals"
    label 'gatk'
    label 'large_mem'
    
    input: 
    path intervals
    path annotated
    path counts
    
    output: 
    path "cohort.filtered.interval_list", emit: filtered_intervals
    
    script:
    """
    echo "${counts.join('\n')}" > counts_files.list
    gatk --java-options "-Xmx16g" FilterIntervals \\
        -L $intervals \\
        --annotated-intervals $annotated \\
        -I counts_files.list \\
        -imr OVERLAPPING_ONLY \\
        -O cohort.filtered.interval_list
    """
}

// Process to determine ploidy for the cohort
process DETERMINE_PLOIDY_COHORT {
    tag "Determine_Ploidy"
    label 'gatk'
    label 'large_mem'
    
    input: 
    path intervals
    path counts
    path priors
    
    output: 
    path "ploidy-calls", emit: ploidy_calls
    
    script:
    """
    echo "${counts.join('\n')}" > counts_files.list
    gatk --java-options "-Xmx16g" DetermineGermlineContigPloidy \\
        -L $intervals \\
        -I counts_files.list \\
        --contig-ploidy-priors $priors \\
        --output . \\
        --output-prefix ploidy \\
        -imr OVERLAPPING_ONLY
    """
}

// Process to scatter intervals
process SCATTER_INTERVALS {
    tag "Scatter_Intervals"
    label 'gatk'
    
    input: 
    path intervals
    
    output: 
    path "temp_*/scattered.interval_list", emit: shards
    
    script:
    """
    gatk --java-options "-Xmx4g" IntervalListTools \\
        --INPUT $intervals \\
        --SUBDIVISION_MODE INTERVAL_COUNT \\
        --SCATTER_CONTENT ${params.scatter_count} \\
        --OUTPUT .
    """
}

// Process germline CNV calling for the cohort
process GERMLINE_CNV_CALLER_COHORT {
    tag "${interval_shard.parent.name}"
    label 'gatk'
    label 'gpu_or_high_cpu'
    
    input: 
    path interval_shard
    path annotated
    path counts
    path ploidy_calls
    
    output: 
    path "shard-calls", emit: calls
    path "shard-model", emit: model
    
    script:
    def sensitivity_params = "--class-coherence-length 1000.0 --cnv-coherence-length 1000.0 --enable-bias-factors false --interval-psi-scale 1.0E-6 --log-mean-bias-standard-deviation 0.01 --sample-psi-scale 1.0E-6"
    """
    echo "${counts.join('\n')}" > counts_files.list
    gatk --java-options "-Xmx16g" GermlineCNVCaller \\
        --run-mode COHORT \\
        -L $interval_shard \\
        -I counts_files.list \\
        --contig-ploidy-calls $ploidy_calls \\
        --annotated-intervals $annotated \\
        --output . \\
        --output-prefix shard \\
        -imr OVERLAPPING_ONLY \\
        $sensitivity_params
    """
}

// Process to post-process CNV calls for individual samples
process POSTPROCESS_CALLS {
    tag "$sample_id"
    label 'gatk'
    publishDir "${outdir}/out_GCNV/vcfs", mode: 'copy', overwrite: true
    
    input: 
    tuple val(index), val(sample_id)
    path model_shards
    path call_shards
    path ploidy_calls
    path dict
    
    output: 
    path "${sample_id}.genotyped_segments.vcf.gz", emit: final_vcf
    path "${sample_id}.genotyped_intervals.vcf.gz", emit: intervals_vcf
    
    script:
    def model_args = model_shards.collect { "--model-shard-path \$it" }.join(" ")
    def call_args = call_shards.collect { "--calls-shard-path \$it" }.join(" ")
    """
    gatk --java-options "-Xmx8g" PostprocessGermlineCNVCalls \\
        \$call_args \\
        \$model_args \\
        --contig-ploidy-calls $ploidy_calls \\
        --sample-index $index \\
        --allosomal-contig chrX --allosomal-contig chrY \\
        --sequence-dictionary $dict \\
        --output-genotyped-intervals ${sample_id}.genotyped_intervals.vcf.gz \\
        --output-genotyped-segments ${sample_id}.genotyped_segments.vcf.gz
    """
}

// Process to compress, sort, index, and annotate each gCNV VCF with TOOL=GATK-gCNV
process BGZIP_SORT_INDEX_VCF {
    tag "${vcf_file.simpleName}"
    label 'gatk'
    publishDir "${outdir}/out_GCNV/vcfs", mode: 'copy', overwrite: true

    input:
    path vcf_file

    output:
    path("*.sorted.vcf.gz"),     emit: sorted_vcf
    path("*.sorted.vcf.gz.tbi"), emit: sorted_vcf_index

    script:
    def sample_name = vcf_file.name.replaceAll(/\.vcf(\.gz)?$/, '')
    def sorted_gz   = "${sample_name}.sorted.vcf.gz"
    """
    # Create extra header line for the TOOL INFO field
    printf '##INFO=<ID=TOOL,Number=1,Type=String,Description="Calling tool">\\n' > extra_header.txt

    # Build a BED annotation file with TOOL=GATK-gCNV for every variant
    bcftools query -f '%CHROM\\t%POS0\\t%END\\n' ${vcf_file} | \\
        awk 'BEGIN{OFS="\\t"} {print \$1, \$2, \$3, "GATK-gCNV"}' | \\
        bgzip -c > ${sample_name}_tool_annot.bed.gz
    tabix -p bed ${sample_name}_tool_annot.bed.gz

    # Annotate the VCF with TOOL=GATK-gCNV in the INFO field
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

workflow GATK_GCNV {
    take:
    bam_ch      // channel: tuple(sample_id, reads, index)
    fasta_ch    // path: reference FASTA
    fai_ch      // path: FASTA index (.fai)
    dict_ch     // path: sequence dictionary (.dict)
    targets_ch  // path: target intervals BED/interval_list

    main:
    // Step 1: Generate ploidy priors from the reference index
    GENERATE_PLOIDY_PRIORS(fai_ch)

    // Step 2: Preprocess intervals with binning and padding
    PREPROCESS_INTERVALS(fasta_ch, fai_ch, dict_ch, targets_ch)

    // Step 3: Annotate intervals with GC content and mappability
    ANNOTATE_INTERVALS(
        PREPROCESS_INTERVALS.out.interval_list,
        fasta_ch, fai_ch, dict_ch
    )

    // Step 4: Collect read counts per sample
    COLLECT_READ_COUNTS(
        bam_ch,
        PREPROCESS_INTERVALS.out.interval_list,
        fasta_ch, fai_ch, dict_ch
    )

    // Step 5: Filter intervals based on annotation and read count statistics
    FILTER_INTERVALS(
        PREPROCESS_INTERVALS.out.interval_list,
        ANNOTATE_INTERVALS.out.annotated_intervals,
        COLLECT_READ_COUNTS.out.counts.map { sample_id, hdf5 -> hdf5 }.collect()
    )

    // Step 6: Determine ploidy for the cohort
    DETERMINE_PLOIDY_COHORT(
        FILTER_INTERVALS.out.filtered_intervals,
        COLLECT_READ_COUNTS.out.counts.map { sample_id, hdf5 -> hdf5 }.collect(),
        GENERATE_PLOIDY_PRIORS.out.priors
    )

    // Step 7: Scatter filtered intervals into shards for parallel processing
    SCATTER_INTERVALS(FILTER_INTERVALS.out.filtered_intervals)

    shards_ch = SCATTER_INTERVALS.out.shards.flatten()

    // Step 8: Run germline CNV caller in cohort mode for each interval shard
    GERMLINE_CNV_CALLER_COHORT(
        shards_ch,
        ANNOTATE_INTERVALS.out.annotated_intervals,
        COLLECT_READ_COUNTS.out.counts.map { sample_id, hdf5 -> hdf5 }.collect(),
        DETERMINE_PLOIDY_COHORT.out.ploidy_calls
    )

    // Step 9: Post-process calls into per-sample genotyped VCFs
    indexed_samples_ch = bam_ch
        .map { sample_id, reads, index -> sample_id }
        .toList()
        .flatMap { samples ->
            samples.withIndex().collect { sample, idx -> tuple(idx, sample) }
        }

    POSTPROCESS_CALLS(
        indexed_samples_ch,
        GERMLINE_CNV_CALLER_COHORT.out.model.collect(),
        GERMLINE_CNV_CALLER_COHORT.out.calls.collect(),
        DETERMINE_PLOIDY_COHORT.out.ploidy_calls,
        dict_ch
    )

    // Step 10: Compress, sort, index, and annotate each segments VCF with TOOL=GATK-gCNV
    BGZIP_SORT_INDEX_VCF(POSTPROCESS_CALLS.out.final_vcf)

    emit:
    sorted_vcf       = BGZIP_SORT_INDEX_VCF.out.sorted_vcf
    sorted_vcf_index = BGZIP_SORT_INDEX_VCF.out.sorted_vcf_index
    intervals_vcf    = POSTPROCESS_CALLS.out.intervals_vcf
}
