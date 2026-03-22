# DRAGEN ICAv2 module (`--workflow dragen`)

## Methodology
This module orchestrates ICAv2-based DRAGEN Germline Enrichment analysis: upload inputs, launch cloud analysis, monitor status, download outputs, then post-process VCFs.

## Processes in this module
1. **UPLOAD_CRAM_FILES** – upload sample CRAM/CRAI (or BAM/BAI) files to ICAv2.
2. **GET_STATIC_FILES** – append reference/bed/static ICA data IDs.
3. **CHECK_FILE_STATUS** – wait until ICA files are ready.
4. **START_ANALYSIS_BATCH** – start DRAGEN batch analysis.
5. **CHECK_ANALYSIS_STATUS** – poll analysis state to completion.
6. **DOWNLOAD_ANALYSIS_OUTPUT** – download output artifacts.
7. **DELETE_DATA** – cleanup ICA-uploaded files/output directory.
8. **ADD_DRAGEN_TOOL_ANNOTATION** – add `TOOL=DRAGEN` annotation to VCFs.
9. **BGZIP_SORT_INDEX_VCF** – sort, compress, and index VCF.
10. **NORMALISE_CNV_QUALITY_SCORES** – normalise caller quality scores.

## Required parameters
- `workflow`: `dragen`
- `projectId`
- `pipelineId`
- `userReference`
- `storageSize`
- `referenceAnalysisDataCode`
- `targetBedAnalysisDataCode`
- `cramAnalysisDataCode`
- `cramIndexAnalysisDataCode`
- `cramReferenceAnalysisDataCode`
- `referenceFileId`
- `cramReferenceFileId`
- `targetBedFileId`
- `localDownloadPath`
- `cramFilePairsUploadPath`
- `icaUploadPath`

## Commonly used optional parameters
- `maxUploadForks`
- `uploadRetries`
- `fileStatusCheckInterval`
- `fileStatusCheckLimit`
- `analysisStatusCheckInterval`
- `analysisStatusCheckLimit`
- `truth_bed`, `probes_bed` (only when running automatic evaluation in `main.nf`)

## Params template
- `params/params-icav2-dragen.json`
