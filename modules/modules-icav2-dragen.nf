#!/usr/bin/env nextflow
nextflow.enable.dsl=2

// =====================================================================================
// PROCESSES FOR ICAv2 DRAGEN GERMLINE ENRICHMENT
// =====================================================================================

outdir = file(params.outdir, type: 'dir')

// Inspired by a workflow by Reagan Cannell @ https://github.com/SBIMB/ica-elwazi/tree/main/nextflow_workflows/cram_input_dragen_ica_workflow

process UPLOAD_CRAM_FILES {
    debug true
    tag "$sampleId"
    label 'icav2-dragen'
    maxForks params.maxUploadForks
    cpus 1
    time '6h'

    errorStrategy { task.exitStatus == 100 ? 'terminate' : 'retry' }
    maxRetries params.uploadRetries

    input:
    tuple val(sampleId), file(cramPair)

    output:
    path "data_upload.txt", emit: dataFile

    script:
    def projectId = params.projectId
    def uploadPath = params.icaUploadPath
    def cramCode = params.cramAnalysisDataCode
    def craiCode = params.cramIndexAnalysisDataCode
    // Support both CRAM/CRAI and BAM/BAI pairs; pair[0] is the primary file, pair[1] the index
    def (cram_file, crai_file) = cramPair
    // Exponential backoff: 15s * 2^(attempt-1) capped at 300s
    def backoffSecs = Math.min(15 * (int)Math.pow(2, task.attempt - 1), 300)
    """
    #!/bin/bash
    set -euo pipefail

    # FORCE CACHE INVALIDATION: Ensure we re-check ICA even if -resume is used
    echo "Verifying upload status for ${sampleId} (attempt ${task.attempt})..."

    # --- JITTER: stagger parallel uploads to avoid thundering-herd on the API ---
    sleep \$((RANDOM % 30))

    time_stamp=\$(date +"%Y-%m-%d %H:%M:%S")
    touch data_upload.txt

    clean_json() {
        sed -n '/^{/,\$p'
    }

    get_or_upload_file() {
        local fpath=\$1
        local folder=\$2
        local pid=\$3
        local fname
        fname=\$(basename "\$fpath")

        # 1. Check existence
        local list_raw
        list_raw=\$(icav2 projectdata list --project-id "\$pid" --parent-folder "\$folder" --file-name "\$fname" -o json 2>/dev/null || true)

        # FIX: Removed "401" from grep to avoid false positives on IDs/Sizes
        if echo "\$list_raw" | grep -qE "Unauthorized|ICA_SEC_002"; then
            echo "CRITICAL ERROR: ICA Authentication Failed during LIST. Check API Key." >&2
            echo "Debug Response: \$list_raw" >&2
            exit 100
        fi

        local existing_id=""
        local existing_status=""
        if [ -n "\$list_raw" ]; then
            existing_id=\$(echo "\$list_raw" | clean_json | jq -r '.items[0].id // empty' 2>/dev/null || true)
            existing_status=\$(echo "\$list_raw" | clean_json | jq -r '.items[0].status // empty' 2>/dev/null || true)
        fi

        if [ -n "\$existing_id" ]; then
            if [ "\$existing_status" == "PARTIAL" ]; then
                echo "WARN: File \$fname exists but is PARTIAL. Deleting and re-uploading..." >&2
                icav2 projectdata delete "\$existing_id" >/dev/null 2>&1 || true
            else
                echo "SKIP: \$fname already exists (ID: \$existing_id)" >&2
                echo "\$existing_id" | tr -d '[:space:]'
                return 0
            fi
        fi
        
        # 2. Upload with extended timeout for large files
        echo "UPLOAD: Uploading \$fname..." >&2
        local upload_raw
        upload_raw=\$(timeout 7200 icav2 projectdata upload "\$fpath" "\$folder" --project-id "\$pid" -o json 2>&1 || true)

        if echo "\$upload_raw" | grep -qE "Unauthorized|ICA_SEC_002"; then
             echo "CRITICAL ERROR: ICA Authentication Failed during UPLOAD. Check API Key." >&2
             echo "Raw Output: \$upload_raw" >&2
             exit 100
        fi
        
        if echo "\$upload_raw" | grep -qE "429|Too Many Requests"; then
             echo "WARNING: Rate limit hit for \$fname. Backing off ${backoffSecs}s before retry." >&2
             sleep ${backoffSecs}
             return 1
        fi

        # Handle transient network/connection errors with backoff
        if echo "\$upload_raw" | grep -qE "connection reset|timeout|Connection refused|Network unreachable|curl: \\("; then
             echo "WARNING: Network error uploading \$fname. Backing off ${backoffSecs}s before retry." >&2
             sleep ${backoffSecs}
             return 1
        fi

        local new_id=""
        new_id=\$(echo "\$upload_raw" | clean_json | jq -r '.id // empty' 2>/dev/null || true)

        if [ -z "\$new_id" ]; then
            new_id=\$(echo "\$upload_raw" | grep -o 'fil.[a-zA-Z0-9]*' | head -1)
        fi

        if [ -z "\$new_id" ]; then
            echo "ERROR: Failed to get ID for uploaded file \$fname" >&2
            echo "DEBUG RAW OUTPUT: \$upload_raw" >&2
            return 1
        fi

        echo "\$new_id" | tr -d '[:space:]'
    }

    echo "[\${time_stamp}] Starting upload process for ${sampleId}..."

    # SEQUENTIAL UPLOAD
    cram_id=\$(get_or_upload_file "${cram_file}" "${uploadPath}" "${projectId}")
    echo "CRAM:\$cram_id" > cram_id.txt

    crai_id=\$(get_or_upload_file "${crai_file}" "${uploadPath}" "${projectId}")
    echo "CRAI:\$crai_id" > crai_id.txt

    cram_final_id=\$(grep 'CRAM:' cram_id.txt | cut -d: -f2 | tr -d '[:space:]')
    crai_final_id=\$(grep 'CRAI:' crai_id.txt | cut -d: -f2 | tr -d '[:space:]')

    if [ -z "\$cram_final_id" ]; then echo "Error: CRAM ID missing"; exit 1; fi
    if [ -z "\$crai_final_id" ]; then echo "Error: CRAI ID missing"; exit 1; fi

    echo "sampleId:${sampleId}" >> data_upload.txt
    echo "${cramCode}:\${cram_final_id}" >> data_upload.txt
    echo "${craiCode}:\${crai_final_id}" >> data_upload.txt

    echo "Processing complete for ${sampleId}."
    """
}

process GET_STATIC_FILES {
    debug true
    tag "Validating Static Files"
    label 'icav2-dragen'
    cpus 1

    input:
    path(dataFile)

    output:
    path "data_static.txt", emit: dataFile

    script:
    def projectId = params.projectId
    def refCode = params.referenceAnalysisDataCode
    def refId = params.referenceFileId
    def bedCode = params.targetBedAnalysisDataCode
    def bedId = params.targetBedFileId
    def cramRefCode = params.cramReferenceAnalysisDataCode
    def cramRefId = params.cramReferenceFileId
    """
    #!/bin/bash
    set -euo pipefail

    cat ${dataFile} > data_static.txt
    echo "" >> data_static.txt
    
    echo "${refCode}:${refId}" >> data_static.txt
    echo "${bedCode}:${bedId}" >> data_static.txt
    echo "${cramRefCode}:${cramRefId}" >> data_static.txt
    
    echo "Static file IDs appended."
    """
}

process CHECK_FILE_STATUS {
    debug true
    tag "Checking Status"
    label 'icav2-dragen'
    cpus 1
    errorStrategy 'retry'
    maxRetries 5
    
    input:
    path(dataFile)

    output:
    path "data_checked.txt", emit: dataFile

    script:
    def interval = params.fileStatusCheckInterval
    def limit = params.fileStatusCheckLimit
    def cramCode = params.cramAnalysisDataCode
    def craiCode = params.cramIndexAnalysisDataCode
    """
    #!/bin/bash
    set -euo pipefail

    echo "DEBUG: Running hardened file status check (large-scale, attempt ${task.attempt})..."
    
    cat ${dataFile} > data_checked.txt

    cram_count=\$(grep -c '^${cramCode}:' ${dataFile} || true)
    crai_count=\$(grep -c '^${craiCode}:' ${dataFile} || true)
    
    echo "Found \$cram_count CRAM/BAM files and \$crai_count CRAI/BAI files to check"

    if [ \$cram_count -eq 0 ]; then
        echo "Error: No primary files (CRAM/BAM) found in input file." >&2
        exit 1
    fi

    if [ \$crai_count -eq 0 ]; then
        echo "Error: No index files (CRAI/BAI) found in input file." >&2
        exit 1
    fi

    mapfile -t cram_ids < <(grep '^${cramCode}:' ${dataFile} | cut -f2- -d: | tr -d '\r')
    mapfile -t crai_ids < <(grep '^${craiCode}:' ${dataFile} | cut -f2- -d: | tr -d '\r')
    
    check_count=0
    cram_auth_failures=0
    crai_auth_failures=0
    
    while true; do
        ((check_count+=1))
        
        echo "[\$(date)]: Attempt \${check_count}/${limit}..."
        
        all_available=true
        
        # Check CRAMs/BAMs
        for i in \${!cram_ids[@]}; do
            cram_id=\${cram_ids[\$i]}
            cram_id=\$(echo "\$cram_id" | tr -d '[:space:]')
            if [ -z "\$cram_id" ]; then continue; fi
            
            cram_status_raw=\$(icav2 projectdata get \${cram_id} -o json 2>&1 || true)
            
            if echo "\$cram_status_raw" | grep -qE "Unauthorized|ICA_SEC_002"; then
                ((cram_auth_failures+=1))
                # Exponential backoff capped at 600s
                raw_backoff=\$(( ${interval} * (1 << (cram_auth_failures - 1)) ))
                backoff=\$(( raw_backoff > 600 ? 600 : raw_backoff ))
                echo "WARNING: Auth/Network glitch checking primary file status (failure #\${cram_auth_failures}). Waiting \${backoff}s..." >&2
                sleep \$backoff
                all_available=false
                continue
            fi
            cram_auth_failures=0
            
            if echo "\$cram_status_raw" | grep -qE "NotFound|not found|404"; then
                 echo "CRITICAL ERROR: File ID \${cram_id} not found in ICA. Please re-run pipeline to re-upload."
                 exit 1
            fi
            
            cram_status=\$(echo "\$cram_status_raw" | sed -n '/^{/,\$p' | jq -r ".details.status" 2>/dev/null || echo "UNKNOWN")
            if [ "\$cram_status" != "AVAILABLE" ]; then
                all_available=false
            fi
            
            # Mini sleep to prevent bursting API; scale down for large batches
            sleep 0.2
        done
        
        # Check CRAIs/BAIs
        for i in \${!crai_ids[@]}; do
            crai_id=\${crai_ids[\$i]}
            crai_id=\$(echo "\$crai_id" | tr -d '[:space:]')
            if [ -z "\$crai_id" ]; then continue; fi

            crai_status_raw=\$(icav2 projectdata get \${crai_id} -o json 2>&1 || true)
            
            if echo "\$crai_status_raw" | grep -qE "Unauthorized|ICA_SEC_002"; then
                ((crai_auth_failures+=1))
                # Exponential backoff capped at 600s
                raw_backoff=\$(( ${interval} * (1 << (crai_auth_failures - 1)) ))
                backoff=\$(( raw_backoff > 600 ? 600 : raw_backoff ))
                echo "WARNING: Auth/Network glitch checking index file status (failure #\${crai_auth_failures}). Waiting \${backoff}s..." >&2
                sleep \$backoff
                all_available=false
                continue
            fi
            crai_auth_failures=0
            
            if echo "\$crai_status_raw" | grep -qE "NotFound|not found|404"; then
                 echo "CRITICAL ERROR: Index ID \${crai_id} not found in ICA. Please re-run pipeline to re-upload."
                 exit 1
            fi
            
            crai_status=\$(echo "\$crai_status_raw" | sed -n '/^{/,\$p' | jq -r ".details.status" 2>/dev/null || echo "UNKNOWN")
            if [ "\$crai_status" != "AVAILABLE" ]; then
                all_available=false
            fi
            
            # Mini sleep to prevent bursting API
            sleep 0.2
        done
        
        if [ "\$all_available" = true ]; then
            echo "✓ All files are AVAILABLE."
            break
        fi
        
        if [ \${check_count} -gt ${limit} ]; then
            echo "✗ Timeout waiting for files."
            exit 1
        fi
        
        sleep ${interval}
    done
    """
}

process START_ANALYSIS_BATCH {
    debug true
    tag "Launch Batch Analysis"
    label 'icav2-dragen'
    cpus 1

    input:
    path(dataFile)

    output:
    path "data_analysis.txt", emit: dataFile

    script:
    def projectId = params.projectId
    def pipelineId = params.pipelineId
    def storageSize = params.storageSize
    def userRefPrefix = params.userReference
    
    def cramCode = params.cramAnalysisDataCode
    def craiCode = params.cramIndexAnalysisDataCode
    def refCode = params.referenceAnalysisDataCode
    def bedCode = params.targetBedAnalysisDataCode
    def cramRefCode = params.cramReferenceAnalysisDataCode
    
    """
    #!/bin/bash
    set -euo pipefail

    cat ${dataFile} > data_analysis.txt
    
    ref_id=\$(grep '^${refCode}:' ${dataFile} | cut -f2- -d: | tr -d '[:space:]')
    bed_id=\$(grep '^${bedCode}:' ${dataFile} | cut -f2- -d: | tr -d '[:space:]')
    cram_ref_id=\$(grep '^${cramRefCode}:' ${dataFile} | cut -f2- -d: | tr -d '[:space:]')

    # Validate minimum 5 samples required for the in-run panel of normals (cnv_enable_in_run_pon)
    sample_count=\$(grep -c "^${cramCode}:" ${dataFile} || true)
    if [ "\$sample_count" -lt 5 ]; then
        echo "Error: DRAGEN Germline Enrichment requires a minimum of 5 samples to enable the in-run panel of normals (cnv_enable_in_run_pon). Found \${sample_count} sample(s); please provide at least 5." >&2
        exit 1
    fi
    echo "Sample count: \${sample_count} (minimum 5 satisfied for in-run panel of normals)"

    # Aggregate lists for batch analysis
    cram_ids_list=\$(grep '^${cramCode}:' ${dataFile} | cut -f2- -d: | tr -d '\r' | sort -u | paste -sd "," -)
    crai_ids_list=\$(grep '^${craiCode}:' ${dataFile} | cut -f2- -d: | tr -d '\r' | sort -u | paste -sd "," -)

    if [ -z "\$cram_ids_list" ]; then echo "Error: No CRAM IDs found"; exit 1; fi
    if [ -z "\$crai_ids_list" ]; then echo "Error: No CRAI IDs found"; exit 1; fi

    echo "Prepared CRAM list: \${cram_ids_list:0:50}..."

    user_ref="${userRefPrefix}-batch"

    echo "Starting batch Nextflow analysis..."

    analysis_response_raw=\$(icav2 projectpipelines start nextflow ${pipelineId} \
        --project-id ${projectId} \
        --user-reference "\${user_ref}" \
        --storage-size ${storageSize} \
        --input ${refCode}:\${ref_id} \
        --input ${bedCode}:\${bed_id} \
        --input ${cramRefCode}:\${cram_ref_id} \
        --input ${cramCode}:\${cram_ids_list} \
        --input ${craiCode}:\${crai_ids_list} \
        --parameters enable_map_align:true \
        --parameters enable_map_align_output:true \
        --parameters output_format:CRAMv3.1 \
        --parameters enable_duplicate_marking:true \
        --parameters enable_variant_caller:true \
        --parameters vc_emit_ref_confidence:GVCF \
        --parameters vc_enable_vcf_output:true \
        --parameters enable_cnv:true \
        --parameters cnv_enable_in_run_pon:true \
        --parameters cnv_enable_gcbias_correction:true \
        --parameters enable_sv:true \
        --parameters qc_coverage_report:cov_report \
        --parameters qc_coverage_filters:"'mapq<20,bq<20'" \
        --parameters qc_detect_contamination:true \
        --parameters enable_variant_annotation:true \
        --parameters samples_per_node:5 2>&1 || true)

    if echo "\$analysis_response_raw" | grep -qE "Unauthorized|ICA_SEC_002"; then
         echo "CRITICAL ERROR: Auth failed launching analysis."
         echo "Raw: \$analysis_response_raw"
         exit 100
    fi

    analysis_id=\$(echo "\$analysis_response_raw" | sed -n '/^{/,\$p' | jq -r ".id" 2>/dev/null || echo "")

    if [ -z "\${analysis_id}" ] || [ "\${analysis_id}" == "null" ]; then
        echo "Error: Failed to launch batch analysis."
        echo "Full API Response: \${analysis_response_raw}"
        exit 1
    fi

    echo "analysisId:\${analysis_id}" >> data_analysis.txt
    echo "analysisRef:\${user_ref}" >> data_analysis.txt
    echo "batchMode:true" >> data_analysis.txt
    """
}

process CHECK_ANALYSIS_STATUS {
    debug true
    tag "Monitor Analysis"
    label 'icav2-dragen'
    cpus 1
    errorStrategy 'retry'
    maxRetries 3

    input:
    path(dataFile)

    output:
    path "data_status.txt", emit: dataFile

    script:
    def interval = params.analysisStatusCheckInterval
    def limit = params.analysisStatusCheckLimit
    """
    #!/bin/bash
    set -euo pipefail

    cat ${dataFile} > data_status.txt

    analysis_id=\$(grep '^analysisId:' ${dataFile} | cut -f2- -d: | tr -d '[:space:]')
    check_count=0
    api_failures=0

    echo "Monitoring analysis \${analysis_id}..."

    while true; do
        ((check_count+=1))

        status_json=\$(icav2 projectanalyses get \${analysis_id} -o json 2>&1 || true)

        # Handle auth / connection failures with exponential backoff
        if echo "\$status_json" | grep -qE "Unauthorized|ICA_SEC_002|connection reset|timeout|curl: \\("; then
            ((api_failures+=1))
            backoff=\$((${interval} * api_failures))
            echo "WARNING: API/network glitch on status check (failure #\${api_failures}). Waiting \${backoff}s..." >&2
            sleep \$backoff
            continue
        fi

        api_failures=0   # reset on successful API call

        status=\$(echo "\$status_json" | sed -n '/^{/,\$p' | jq -r ".status // empty" 2>/dev/null | tr -d '[:space:]')

        if [ -z "\$status" ]; then
            echo "WARNING: Could not parse status from API response. Retrying in ${interval}s..." >&2
            sleep ${interval}
            continue
        fi

        echo "[\$(date)]: Status is \${status}"

        if [[ "\${status}" == "SUCCEEDED" ]]; then
            echo "analysisStatus:SUCCEEDED" >> data_status.txt
            break
        elif [[ "\${status}" == "FAILED" || "\${status}" == "FAILED_FINAL" || "\${status}" == "ABORTED" ]]; then
            echo "analysisStatus:FAILED" >> data_status.txt
            break
        fi

        if [ \${check_count} -gt ${limit} ]; then
            echo "Timeout waiting for analysis to complete."
            echo "analysisStatus:TIMEOUT" >> data_status.txt
            break
        fi

        sleep ${interval}
    done
    """
}

process DOWNLOAD_ANALYSIS_OUTPUT {
    debug true
    tag "Download Output"
    label 'icav2-dragen'
    cpus 1

    input:
    path(dataFile)

    output:
    path "data_final.txt", emit: dataFile

    script:
    def downloadPath = "${outdir}/out_DRAGEN"
    """
    #!/bin/bash
    set -euo pipefail

    cat ${dataFile} > data_final.txt

    status=\$(grep '^analysisStatus:' ${dataFile} | cut -f2- -d:)
    analysis_id=\$(grep '^analysisId:' ${dataFile} | cut -f2- -d:)

    if [ "\$status" == "SUCCEEDED" ]; then
        echo "Fetching output folder info..."
        output_json=\$(icav2 projectanalyses output \${analysis_id} -o json)
        folder_id=\$(echo \${output_json} | jq -r ".items[0].data[0].dataId")

        if [ -n "\${folder_id}" ] && [ "\${folder_id}" != "null" ]; then
            mkdir -p "${downloadPath}"
            echo "Downloading output folder \${folder_id} to ${downloadPath}..."
            icav2 projectdata download \${folder_id} "${downloadPath}"
            echo "outputFolderId:\${folder_id}" >> data_final.txt
            echo "deleteData:true" >> data_final.txt
        else
            echo "Warning: No output folder ID found in analysis response."
            echo "deleteData:false" >> data_final.txt
        fi
    else
        echo "Analysis failed or timed out. Skipping download."
        echo "deleteData:false" >> data_final.txt
    fi
    """
}

process DELETE_DATA {
    debug true
    tag "Cleanup"
    label 'icav2-dragen'
    cpus 1
    errorStrategy 'ignore'

    input:
    path(dataFile)

    output:
    stdout

    script:
    def cramCode = params.cramAnalysisDataCode
    def craiCode = params.cramIndexAnalysisDataCode
    """
    #!/bin/bash

    delete_flag=\$(grep '^deleteData:' ${dataFile} | cut -f2- -d: | tr -d '[:space:]')

    if [ "\$delete_flag" == "true" ]; then
        # Delete ALL uploaded primary files (CRAM or BAM) for every sample
        while IFS= read -r line; do
            file_id=\$(echo "\$line" | cut -d: -f2 | tr -d '[:space:]')
            if [ -n "\$file_id" ]; then
                echo "Deleting primary file \${file_id}..."
                icav2 projectdata delete "\${file_id}" || true
            fi
        done < <(grep "^${cramCode}:" "${dataFile}")

        # Delete ALL uploaded index files (CRAI or BAI) for every sample
        while IFS= read -r line; do
            idx_id=\$(echo "\$line" | cut -d: -f2 | tr -d '[:space:]')
            if [ -n "\$idx_id" ]; then
                echo "Deleting index file \${idx_id}..."
                icav2 projectdata delete "\${idx_id}" || true
            fi
        done < <(grep "^${craiCode}:" "${dataFile}")

        folder_id=\$(grep '^outputFolderId:' ${dataFile} | cut -f2- -d: | tr -d '[:space:]')
        if [ -n "\${folder_id}" ]; then
            echo "Deleting analysis output folder \${folder_id}..."
            icav2 projectdata delete "\${folder_id}" || true
        fi

        echo "Cleanup complete."
    else
        echo "Skipping cleanup (Analysis failed or download skipped)."
    fi
    """
}

process ADD_DRAGEN_TOOL_ANNOTATION {
    debug true
    tag "Add TOOL=DRAGEN annotation"
    label 'bcftools'
    publishDir "${outdir}/out_DRAGEN/annotated_vcfs", mode: 'copy', overwrite: true
    cpus 1

    input:
    path(dataFile)

    output:
    path "*_DRAGEN.annotated.vcf.gz", emit: annotated_vcfs, optional: true

    script:
    def downloadPath = "${outdir}/out_DRAGEN"
    """
    #!/bin/bash
    set -euo pipefail

    echo "Annotating downloaded DRAGEN VCF files with TOOL=DRAGEN in INFO field..."

    printf '##INFO=<ID=TOOL,Number=1,Type=String,Description="Calling tool">\\n' > extra_header.txt

    found=0
    while IFS= read -r vcf; do
        found=1
        sample_name=\$(basename "\$vcf" | sed 's/\\.cnv\\.vcf\\.gz\$//' | sed 's/\\.vcf\\.gz\$//' | sed 's/\\.vcf\$//')

        bcftools query -f '%CHROM\\t%POS0\\t%END\\n' "\$vcf" | \\
            awk 'BEGIN{OFS="\\t"} {print \$1, \$2, \$3, "DRAGEN"}' | \\
            bgzip -c > "\${sample_name}_tool_annot.bed.gz"
        tabix -p bed "\${sample_name}_tool_annot.bed.gz"

        bcftools annotate \\
            -a "\${sample_name}_tool_annot.bed.gz" \\
            -c CHROM,FROM,TO,INFO/TOOL \\
            -h extra_header.txt \\
            "\$vcf" \\
            -O z -o "\${sample_name}_DRAGEN.annotated.vcf.gz"

        rm -f "\${sample_name}_tool_annot.bed.gz" "\${sample_name}_tool_annot.bed.gz.tbi"
    done < <(find ${downloadPath} -type f -name "*.cnv.vcf.gz" 2>/dev/null || true)

    if [ \$found -eq 0 ]; then
        echo "No VCF files found in ${downloadPath}. Skipping annotation."
    else
        echo "TOOL=DRAGEN annotation complete."
    fi
    """
}

// Process to compress, sort, and index each annotated DRAGEN VCF
process BGZIP_SORT_INDEX_VCF {
    tag "${vcf_file.simpleName}"
    label 'bcftools'
    publishDir "${outdir}/out_DRAGEN/sorted_vcfs", mode: 'copy', overwrite: true

    input:
    path vcf_file

    output:
    path("*.sorted.vcf.gz"),     emit: sorted_vcf
    path("*.sorted.vcf.gz.tbi"), emit: sorted_vcf_index

    script:
    def sample_name = vcf_file.name - '.annotated.vcf.gz'
    def sorted_gz   = "${sample_name}.sorted.vcf.gz"
    """
    bcftools sort ${vcf_file} -o ${sorted_gz} -O z
    tabix -p vcf ${sorted_gz}
    """
}

// Process to normalise CNV quality scores to a common scale
process NORMALISE_CNV_QUALITY_SCORES {
    tag "${vcf.simpleName}"
    label 'pysam'
    publishDir "${outdir}/out_DRAGEN/vcfs", mode: 'copy', overwrite: true

    input:
    path vcf

    output:
    path("*.normalised.vcf.gz"),     emit: normalised_vcf
    path("*.normalised.vcf.gz.tbi"), emit: normalised_vcf_index

    script:
    def sample_name = vcf.name - '.sorted.vcf.gz'
    def normalised_gz = "${sample_name}.normalised.vcf.gz"
    """
    normalise_cnv_caller_quality_scores.py \\
        --input_vcf ${vcf} \\
        --output_vcf ${sample_name}.normalised.vcf \\
        --caller DRAGEN
    bgzip -c ${sample_name}.normalised.vcf > ${normalised_gz}
    tabix -p vcf ${normalised_gz}
    """
}

// =====================================================================================
// SUB-WORKFLOW TO CHAIN THE PROCESSES TOGETHER
// =====================================================================================

workflow DRAGEN {
    take:
    cram_ch  // channel: tuple(sampleId, [cram_file, crai_file])

    main:
    // Validate minimum 5 samples required for the in-run panel of normals (cnv_enable_in_run_pon)
    validated_cram_ch = cram_ch
        .collect()
        .flatMap { items ->
            if (items.size() < 5) {
                error "DRAGEN Germline Enrichment requires a minimum of 5 samples to enable the in-run panel of normals (cnv_enable_in_run_pon). Found ${items.size()} sample(s); please provide at least 5."
            }
            return items
        }

    // Step 1: Upload CRAM and CRAI files to ICA for each sample
    UPLOAD_CRAM_FILES(validated_cram_ch)

    // Step 2: Combine all per-sample upload data into a single file
    combined_upload_ch = UPLOAD_CRAM_FILES.out.dataFile
        .collectFile(name: 'combined_data.txt')

    // Step 3: Append static reference file IDs to the combined data file
    GET_STATIC_FILES(combined_upload_ch)

    // Step 4: Verify all uploaded files are available in ICA
    CHECK_FILE_STATUS(GET_STATIC_FILES.out.dataFile)

    // Step 5: Launch the DRAGEN batch analysis on ICA
    START_ANALYSIS_BATCH(CHECK_FILE_STATUS.out.dataFile)

    // Step 6: Monitor the analysis until completion
    CHECK_ANALYSIS_STATUS(START_ANALYSIS_BATCH.out.dataFile)

    // Step 7: Download the analysis output from ICA
    DOWNLOAD_ANALYSIS_OUTPUT(CHECK_ANALYSIS_STATUS.out.dataFile)

    // Step 8: Clean up temporary CRAM files and output folder from ICA
    DELETE_DATA(DOWNLOAD_ANALYSIS_OUTPUT.out.dataFile)

    // Step 9: Annotate downloaded VCF files with TOOL=DRAGEN in INFO field
    ADD_DRAGEN_TOOL_ANNOTATION(DOWNLOAD_ANALYSIS_OUTPUT.out.dataFile)

    // Step 10: Sort and index each annotated VCF
    BGZIP_SORT_INDEX_VCF(ADD_DRAGEN_TOOL_ANNOTATION.out.annotated_vcfs.flatten())

    // Step 11: Normalise quality scores to a common scale
    NORMALISE_CNV_QUALITY_SCORES(BGZIP_SORT_INDEX_VCF.out.sorted_vcf.flatten())

    emit:
    result               = DOWNLOAD_ANALYSIS_OUTPUT.out.dataFile
    annotated_vcfs       = ADD_DRAGEN_TOOL_ANNOTATION.out.annotated_vcfs
    sorted_vcf           = BGZIP_SORT_INDEX_VCF.out.sorted_vcf
    sorted_vcf_index     = BGZIP_SORT_INDEX_VCF.out.sorted_vcf_index
    normalised_vcf       = NORMALISE_CNV_QUALITY_SCORES.out.normalised_vcf
    normalised_vcf_index = NORMALISE_CNV_QUALITY_SCORES.out.normalised_vcf_index
}
