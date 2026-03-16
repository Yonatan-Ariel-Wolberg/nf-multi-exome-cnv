#!/usr/bin/env nextflow
nextflow.enable.dsl=2

// =============================================================================
// INDELIBLE MODULE
// =============================================================================

// REQUIRED FILES
ref               = file(params.ref, type: 'file')
priors            = file(params.priors, type: 'file')
indelible_conf    = file(params.indelible_conf, type: 'file')
outdir            = file(params.outdir, type: 'dir')

// 1. THE FETCH COMMAND EXTRACTS THE READS FROM THE BAM FILE, IT TAKES 2 ARGUMENTS:
process run_Fetch {
    tag { sample }
    label 'indelible'
    
    input:
    tuple val(sample), path(bam), path(bai)
    
    output:
    tuple val(sample), path("${bam}.sc_reads"), path(bam), path(bai), emit: sc_reads
    
    script:
    """
    export REF_PATH=/dataB/aux/38/samtools_ref_cache/%2s/%2s/%s:http://www.ebi.ac.uk/ena/cram/md5/%s
    export REF_CACHE=/dataB/aux/38/samtools_ref_cache/%2s/%2s/%s
    indelible.py fetch --config ${indelible_conf} --i ${bam} --o ${bam}.sc_reads
    """
}

// 2. THE AGGREGATE MERGES INFORMATION ACROSS READS TOWARDS A POSITION-LEVEL VIEW OF THE DATA:
process run_Aggregate {
    tag { sample }
    label 'indelible'
    
    input:
    tuple val(sample), path(sc_read), path(bam), path(bai)
    
    output:
    tuple val(sample), path("${bam}.counts"), emit: counts
    
    script:
    """
    indelible.py aggregate --i ${sc_read} --b ${bam} --o ${bam}.counts --r ${ref} --config ${indelible_conf}
    """
}

// 3. THE SCORE COMMAND SCORES POSITIONS BASED ON THE READ INFORMATION AND SEQUENCE CONTEXT:
process run_Score {
    tag { sample }
    label 'indelible'
    
    input:
    tuple val(sample), path(count)

    output:
    tuple val(sample), path("${count}.scored"), emit: scores
    path("${count}.scored"), emit: database_in
    
    script:
    """
    indelible.py score --i ${count} --o ${count}.scored --config ${indelible_conf}
    """
}

// 4. THE DATABASE COMMAND GENERATES THE ALLELE FREQUENCY AND BREAKPOINT DATABASE REQUIRED FOR THE NEXT STEP – ANNOTATE
process run_Database {
    tag { "Indel_DB" }
    label 'indelible'
    cpus 6
    publishDir "${outdir}/out_INDELIBLE/database", mode: 'copy', overwrite: true

    input:
    path(score)

    output:
    path("InDelible_db.tsv"), emit: indel_database
    
    script:
    """
    ls *.scored > scores.txt
    indelible.py database --f scores.txt --o InDelible_db.tsv --r ${ref} --priors ${priors} --config ${indelible_conf} --tb ${task.cpus}
    """
}

// 5. THE ANNOTATE COMMAND ENRICHES THE RESULT WITH GENE/EXON ANNOTATIONS AND MERGES THE DATABASE RESULTS WITH THE POSITION FILE:
process run_Annotate {
    tag { "${score}" }
    label 'indelible'
    publishDir "${outdir}/out_INDELIBLE/annotations", mode: 'copy', overwrite: true

    input:
    path(database)
    tuple val(sample), path(score)

    output:
    tuple val(sample), path("${score}.annotated"), emit: annotated
    
    script:
    """
    indelible.py annotate --i ${score} --o ${score}.annotated --d ${database} --config ${indelible_conf}
    """
}

// 6. ONE CAN THEN LOOK FOR DE NOVO MUTATION EVENTS USING THE DENOVO COMMAND:
// TRIO
process run_DenovoTrio {
    tag { sample }
    label 'indelible'
    publishDir "${outdir}/out_INDELIBLE/denovo_annotation_trio", mode: 'copy', overwrite: true

    input:
    tuple val(sample), path(child_bam), path(child_bai), path(mom_bam), path(mom_bai), path(dad_bam), path(dad_bai), path(annotation)
    
    output:
    tuple val(sample), path("${annotation}.denovo.tsv"), emit: indelible_denovo

    script:
    """
    indelible.py denovo --c ${annotation} --m ${mom_bam} --p ${dad_bam} --o ${annotation}.denovo.tsv --config ${indelible_conf}
    """    
}

// MOM
process run_DenovoMom {
    tag { sample }
    label 'indelible'
    publishDir "${outdir}/out_INDELIBLE/denovo_annotation_mom", mode: 'copy', overwrite: true

    input:
    tuple val(sample), path(child_bam), path(child_bai), path(mom_bam), path(mom_bai), path(annotation)
    
    output:
    tuple val(sample), path("${annotation}.denovo.tsv"), emit: indelible_denovo_mom

    script:
    """
    indelible.py denovo --c ${annotation} --m ${mom_bam} --o ${annotation}.denovo.tsv --config ${indelible_conf}
    """    
}

// DAD
process run_DenovoDad {
    label 'indelible'
    tag { sample }

    publishDir "${outdir}/out_INDELIBLE/denovo_annotation_dad", mode: 'copy', overwrite: true

    input:
    tuple val(sample), path(child_bam), path(child_bai), path(dad_bam), path(dad_bai), path(annotation)
    
    output:
    tuple val(sample), path("${annotation}.denovo.tsv"), emit: indelible_denovo_dad

    script:
    """
    indelible.py denovo --c ${annotation} --p ${dad_bam} --o ${annotation}.denovo.tsv --config ${indelible_conf}
    """    
}

// Filter the results from Indelible
process filterINDELIBLE {
    tag { 'filter_cnvs' }
    label 'indelible'
    publishDir "${outdir}/out_INDELIBLE/filtered", mode: 'copy', overwrite: true
    
    input:
    tuple val(sample), path(annotated)
    
    output:
    path("${sample}.annotated.filtered.tsv"), emit: filtered_cnvs
    
    script:
    """
    awk '{ if ((\$39 < 2) && (\$40 < 2)) { print } }' ${annotated} > ${sample}.annotated.filtered.tsv
    """
}

// Process to compress, sort, index, and annotate each INDELIBLE VCF with TOOL=INDELIBLE
process BGZIP_SORT_INDEX_VCF {
    tag "${vcf_file.simpleName}"
    label 'bcftools'
    publishDir "${outdir}/out_INDELIBLE/vcfs", mode: 'copy', overwrite: true

    input:
    path vcf_file

    output:
    path("*.sorted.vcf.gz"),     emit: sorted_vcf
    path("*.sorted.vcf.gz.tbi"), emit: sorted_vcf_index

    script:
    def sample_name = vcf_file.simpleName
    def sorted_gz   = "${sample_name}.sorted.vcf.gz"
    """
    # Create extra header line for the TOOL INFO field
    printf '##INFO=<ID=TOOL,Number=1,Type=String,Description="Calling tool">\\n' > extra_header.txt

    # Build a BED annotation file with TOOL=INDELIBLE for every variant
    bcftools query -f '%CHROM\\t%POS0\\t%END\\n' ${vcf_file} | \\
        awk 'BEGIN{OFS="\\t"} {print \$1, \$2, \$3, "INDELIBLE"}' | \\
        bgzip -c > ${sample_name}_tool_annot.bed.gz
    tabix -p bed ${sample_name}_tool_annot.bed.gz

    # Annotate the VCF with TOOL=INDELIBLE in the INFO field
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

workflow INDELIBLE {
    take:
    crams_ch        // channel: tuple(sample, bam, bai)
    cram_trios_ch   // channel: tuple(sample, child_bam, child_bai, mom_bam, mom_bai, dad_bam, dad_bai)
    cram_mom_ch     // channel: tuple(sample, child_bam, child_bai, mom_bam, mom_bai)
    cram_dad_ch     // channel: tuple(sample, child_bam, child_bai, dad_bam, dad_bai)

    main:
    // Step 1: Extract split/clipped reads from BAM files
    run_Fetch(crams_ch)

    // Step 2: Aggregate read information to a position-level view
    run_Aggregate(run_Fetch.out.sc_reads)

    // Step 3: Score positions based on read information and sequence context
    run_Score(run_Aggregate.out.counts)

    // Step 4: Generate the allele frequency and breakpoint database
    run_Database(run_Score.out.database_in.collect())

    // Step 5: Annotate positions with gene/exon information
    run_Annotate(run_Database.out.indel_database, run_Score.out.scores)

    // Step 6: Identify de novo mutations in complete trios
    run_DenovoTrio(cram_trios_ch.join(run_Annotate.out.annotated))

    // Step 7: Identify de novo mutations with mother only
    run_DenovoMom(cram_mom_ch.join(run_Annotate.out.annotated))

    // Step 8: Identify de novo mutations with father only
    run_DenovoDad(cram_dad_ch.join(run_Annotate.out.annotated))

    // Step 9: Filter annotated results by confidence thresholds
    filterINDELIBLE(run_Annotate.out.annotated)

    emit:
    filtered_cnvs        = filterINDELIBLE.out.filtered_cnvs
    indelible_denovo     = run_DenovoTrio.out.indelible_denovo
    indelible_denovo_mom = run_DenovoMom.out.indelible_denovo_mom
    indelible_denovo_dad = run_DenovoDad.out.indelible_denovo_dad
}
