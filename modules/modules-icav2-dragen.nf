#!/usr/bin/env nextflow
nextflow.enable.dsl=2

// =====================================================================================
// PROCESSES FOR ICAv2 DRAGEN GERMLINE ENRICHMENT
// =====================================================================================

// Inspired by a workflow by Reagan Cannell @ https://github.com/SBIMB/ica-elwazi/tree/main/nextflow_workflows/cram_input_dragen_ica_workflow

process uploadCramFiles {
    debug true
    tag "$sampleId"
    label 'icav2-dragen'
    maxForks params.maxUploadForks
    cpus 1
    time '4h'

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
    def (cram_file, crai_file) = cramPair
    """
    #!/bin/bash
    set -euo pipefail

    # FORCE CACHE INVALIDATION: Ensure we re-check ICA even if -resume is used
    echo "Verifying upload status for ${sampleId} (v6)..."

    # --- JITTER ---
    sleep \$((RANDOM % 15))

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
        
        # 2. Upload
        echo "UPLOAD: Uploading \$fname..." >&2
        local upload_raw
        upload_raw=\$(timeout 3600 icav2 projectdata upload "\$fpath" "\$folder" --project-id "\$pid" -o json 2>&1 || true)

        if echo "\$upload_raw" | grep -qE "Unauthorized|ICA_SEC_002"; then
             echo "CRITICAL ERROR: ICA Authentication Failed during UPLOAD. Check API Key." >&2
             echo "Raw Output: \$upload_raw" >&2
             exit 100
        fi
        
        if echo "\$upload_raw" | grep -qE "429|Too Many Requests"; then
             echo "WARNING: Rate limit hit for \$fname. Failing task to trigger retry." >&2
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

process getStaticFiles {
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

process checkFileStatus {
    debug true
    tag "Checking Status"
    label 'icav2-dragen'
    cpus 1
    errorStrategy 'retry'
    maxRetries 2
    
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

    echo "DEBUG: Running hardened file status check (v6 - fixed ID parsing)..."
    
    cat ${dataFile} > data_checked.txt

    cram_count=\$(grep -c '^${cramCode}:' ${dataFile} || true)
    crai_count=\$(grep -c '^${craiCode}:' ${dataFile} || true)
    
    echo "Found \$cram_count CRAM files and \$crai_count CRAI files to check"

    if [ \$cram_count -eq 0 ]; then
        echo "Error: No CRAM files found in input file." >&2
        exit 1
    fi

    if [ \$crai_count -eq 0 ]; then
        echo "Error: No CRAI files found in input file." >&2
        exit 1
    fi

    mapfile -t cram_ids < <(grep '^${cramCode}:' ${dataFile} | cut -f2- -d: | tr -d '\r')
    mapfile -t crai_ids < <(grep '^${craiCode}:' ${dataFile} | cut -f2- -d: | tr -d '\r')
    
    check_count=0
    
    while true; do
        ((check_count+=1))
        
        echo "[\$(date)]: Attempt \${check_count}/${limit}..."
        
        all_available=true
        
        # Check CRAMs
        for i in \${!cram_ids[@]}; do
            cram_id=\${cram_ids[\$i]}
            cram_id=\$(echo "\$cram_id" | tr -d '[:space:]')
            if [ -z "\$cram_id" ]; then continue; fi
            
            cram_status_raw=\$(icav2 projectdata get \${cram_id} -o json 2>&1 || true)
            
            if echo "\$cram_status_raw" | grep -qE "Unauthorized|ICA_SEC_002"; then
                 echo "WARNING: Auth/Network glitch checking CRAM status. Retrying in ${interval}s..."
                 all_available=false
                 continue
            fi
            
            if echo "\$cram_status_raw" | grep -qE "NotFound|not found|404"; then
                 echo "CRITICAL ERROR: CRAM ID \${cram_id} not found in ICA. Please re-run pipeline to re-upload."
                 exit 1
            fi
            
            cram_status=\$(echo "\$cram_status_raw" | sed -n '/^{/,\$p' | jq -r ".details.status" 2>/dev/null || echo "UNKNOWN")
            if [ "\$cram_status" != "AVAILABLE" ]; then
                all_available=false
            fi
            
            # Mini sleep to prevent bursting API
            sleep 0.5 
        done
        
        # Check CRAIs
        for i in \${!crai_ids[@]}; do
            crai_id=\${crai_ids[\$i]}
            crai_id=\$(echo "\$crai_id" | tr -d '[:space:]')
            if [ -z "\$crai_id" ]; then continue; fi

            crai_status_raw=\$(icav2 projectdata get \${crai_id} -o json 2>&1 || true)
            
            if echo "\$crai_status_raw" | grep -qE "Unauthorized|ICA_SEC_002"; then
                 echo "WARNING: Auth/Network glitch checking CRAI status. Retrying in ${interval}s..."
                 all_available=false
                 continue
            fi
            
            if echo "\$crai_status_raw" | grep -qE "NotFound|not found|404"; then
                 echo "CRITICAL ERROR: CRAI ID \${crai_id} not found in ICA. Please re-run pipeline to re-upload."
                 exit 1
            fi
            
            crai_status=\$(echo "\$crai_status_raw" | sed -n '/^{/,\$p' | jq -r ".details.status" 2>/dev/null || echo "UNKNOWN")
            if [ "\$crai_status" != "AVAILABLE" ]; then
                all_available=false
            fi
            
            # Mini sleep to prevent bursting API
            sleep 0.5
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

process startAnalysisBatch {
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
        --parameters qc_coverage_filters:"'mapq<1'" \
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

process checkAnalysisStatus {
    debug true
    tag "Monitor Analysis"
    label 'icav2-dragen'
    cpus 1

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

    echo "Monitoring analysis \${analysis_id}..."

    while true; do
        ((check_count+=1))

        status_json=\$(icav2 projectanalyses get \${analysis_id} -o json)
        status=\$(echo \${status_json} | jq -r ".status" | tr -d '[:space:]')

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

process downloadAnalysisOutput {
    debug true
    tag "Download Output"
    label 'icav2-dragen'
    cpus 1

    input:
    path(dataFile)

    output:
    path "data_final.txt", emit: dataFile

    script:
    def downloadPath = params.localDownloadPath
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

process deleteData {
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

    delete_flag=\$(grep '^deleteData:' ${dataFile} | cut -f2- -d:)

    if [ "\$delete_flag" == "true" ]; then
        cram_id=\$(grep '^${cramCode}:' ${dataFile} | cut -f2- -d:)
        crai_id=\$(grep '^${craiCode}:' ${dataFile} | cut -f2- -d:)
        folder_id=\$(grep '^outputFolderId:' ${dataFile} | cut -f2- -d:)

        echo "Deleting temporary CRAM file \${cram_id}..."
        icav2 projectdata delete \${cram_id} || true

        echo "Deleting temporary CRAI file \${crai_id}..."
        icav2 projectdata delete \${crai_id} || true

        echo "Deleting analysis output folder \${folder_id}..."
        icav2 projectdata delete \${folder_id} || true

        echo "Cleanup complete."
    else
        echo "Skipping cleanup (Analysis failed or download skipped)."
    fi
    """
}

process addDragenToolAnnotation {
    debug true
    tag "Add TOOLS=DRAGEN annotation"
    label 'bcftools'
    publishDir "${params.localDownloadPath}", mode: 'copy', overwrite: true
    cpus 1

    input:
    path(dataFile)

    output:
    path "*_DRAGEN_output.annotated.vcf.gz", emit: annotated_vcfs, optional: true

    script:
    def downloadPath = params.localDownloadPath
    """
    #!/bin/bash
    set -euo pipefail

    echo "Annotating downloaded DRAGEN VCF files with TOOLS=DRAGEN in INFO field..."

    printf '##INFO=<ID=TOOLS,Number=1,Type=String,Description="Source caller tool">\\n' > extra_header.txt

    found=0
    while IFS= read -r vcf; do
        found=1
        sample_name=\$(basename "\$vcf" | sed 's/\\.vcf\\.gz\$//' | sed 's/\\.vcf\$//')

        bcftools query -f '%CHROM\\t%POS0\\t%END\\n' "\$vcf" | \\
            awk 'BEGIN{OFS="\\t"} {print \$1, \$2, \$3, "DRAGEN"}' | \\
            bgzip -c > "\${sample_name}_tools_annot.bed.gz"
        tabix -p bed "\${sample_name}_tools_annot.bed.gz"

        bcftools annotate \\
            -a "\${sample_name}_tools_annot.bed.gz" \\
            -c CHROM,FROM,TO,INFO/TOOLS \\
            -h extra_header.txt \\
            "\$vcf" \\
            -O z -o "\${sample_name}_DRAGEN_output.annotated.vcf.gz"

        rm -f "\${sample_name}_tools_annot.bed.gz" "\${sample_name}_tools_annot.bed.gz.tbi"
    done < <(find ${downloadPath} -type f \\( -name "*.vcf.gz" -o -name "*.vcf" \\) 2>/dev/null || true)

    if [ \$found -eq 0 ]; then
        echo "No VCF files found in ${downloadPath}. Skipping annotation."
    else
        echo "TOOLS=DRAGEN annotation complete."
    fi
    """
}

// =====================================================================================
// SUB-WORKFLOW TO CHAIN THE PROCESSES TOGETHER
// =====================================================================================

workflow DRAGEN {
    take:
    cram_ch  // channel: tuple(sampleId, [cram_file, crai_file])

    main:
    // Step 1: Upload CRAM and CRAI files to ICA for each sample
    uploadCramFiles(cram_ch)

    // Step 2: Combine all per-sample upload data into a single file
    combined_upload_ch = uploadCramFiles.out.dataFile
        .collectFile(name: 'combined_data.txt')

    // Step 3: Append static reference file IDs to the combined data file
    getStaticFiles(combined_upload_ch)

    // Step 4: Verify all uploaded files are available in ICA
    checkFileStatus(getStaticFiles.out.dataFile)

    // Step 5: Launch the DRAGEN batch analysis on ICA
    startAnalysisBatch(checkFileStatus.out.dataFile)

    // Step 6: Monitor the analysis until completion
    checkAnalysisStatus(startAnalysisBatch.out.dataFile)

    // Step 7: Download the analysis output from ICA
    downloadAnalysisOutput(checkAnalysisStatus.out.dataFile)

    // Step 8: Clean up temporary CRAM files and output folder from ICA
    deleteData(downloadAnalysisOutput.out.dataFile)

    // Step 9: Annotate downloaded VCF files with TOOLS=DRAGEN in INFO field
    addDragenToolAnnotation(downloadAnalysisOutput.out.dataFile)

    emit:
    result = downloadAnalysisOutput.out.dataFile
    annotated_vcfs = addDragenToolAnnotation.out.annotated_vcfs
}
