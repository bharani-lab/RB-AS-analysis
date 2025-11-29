#!/bin/bash
#
# SCRIPT 1: SRA DATA DOWNLOAD AND PREPROCESSING
#
# Purpose: Download RNA-seq FASTQ files from NCBI SRA for retinoblastoma tumors and retinal controls
# Output: FASTQ files organized in data/raw/ directory
# Dependencies: fastq-dump (SRA Toolkit), conda or manual installation
#
# Author: A. Mohamed Hameed Aslam (AMRF Lab, Allagappa University)
# Date: 2025-11-30
#

set -euo pipefail

# ================================
# LOGGING AND UTILITY FUNCTIONS
# ================================

logmessage() {
    local msg="$1"
    local timestamp=$(date '+%Y-%m-%d %H:%M:%S')
    local logentry="[$timestamp] $msg"
    echo "$logentry"
    echo "$logentry" >> "$DOWNLOAD_LOG"
}

# ================================
# DIRECTORY CONFIGURATION
# ================================

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
PROJECT_ROOT="${SCRIPT_DIR%/scripts}"
RAW_DATA_DIR="${PROJECT_ROOT}/data/raw"
LOG_DIR="${PROJECT_ROOT}/logs"
DOWNLOAD_LOG="${LOG_DIR}/download_$(date +%Y%m%d_%H%M%S).txt"

mkdir -p "${RAW_DATA_DIR}/RB"
mkdir -p "${RAW_DATA_DIR}/Control"
mkdir -p "${LOG_DIR}"

logmessage "===================================================="
logmessage "RNA-seq SRA Data Download Pipeline"
logmessage "===================================================="
logmessage "Project Root: $PROJECT_ROOT"
logmessage "Output Directory: $RAW_DATA_DIR"
logmessage "Log File: $DOWNLOAD_LOG"

# ================================
# DEPENDENCY CHECK
# ================================

if ! command -v fastq-dump &> /dev/null; then
    logmessage "ERROR: fastq-dump not found. Please install SRA Toolkit."
    logmessage "Installation: conda install -c bioconda sra-tools"
    exit 1
fi

logmessage "SRA Toolkit found: $(fastq-dump --version 2>&1 | head -1)"

# ================================
# SAMPLE ACCESSIONS
# ================================

# RB TUMOR SAMPLES (SRA accessions)
RB_ACCESSIONS=(
    "SRR13025501"  # RB sample 1
    "SRR13025502"  # RB sample 2
    "SRR13025503"  # RB sample 3
    "SRR13025504"  # RB sample 4
    "SRR13025505"  # RB sample 5
    # Add more RB accessions as needed (production: 50 samples)
)

# CONTROL SAMPLES - Normal retinal tissue
CONTROL_ACCESSIONS=(
    "SRR11292097"  # Control sample 1
    "SRR11292098"  # Control sample 2
    "SRR11292099"  # Control sample 3
    "SRR11292100"  # Control sample 4
    # Add more control accessions as needed (production: 17 samples)
)

logmessage "Configuration:"
logmessage "  - RB samples: ${#RB_ACCESSIONS[@]}"
logmessage "  - Control samples: ${#CONTROL_ACCESSIONS[@]}"

# ================================
# DOWNLOAD FUNCTION
# ================================

download_and_convert_sample() {
    local accession="$1"
    local sample_type="$2"  # "RB" or "Control"
    local output_dir="${RAW_DATA_DIR}/${sample_type}"
    
    mkdir -p "$output_dir"
    
    logmessage "Downloading [$sample_type] $accession..."
    
    if fastq-dump \
        --defline-seq '@$sn[_$rn]/$ri' \
        --defline-qual '+' \
        --split-files \
        --gzip \
        --outdir "$output_dir" \
        --skip-technical \
        "$accession" 2>> "$DOWNLOAD_LOG"; then
        logmessage "✓ Successfully downloaded: $accession"
        return 0
    else
        logmessage "✗ FAILED to download: $accession - Check log for details"
        return 1
    fi
}

# ================================
# MAIN DOWNLOAD LOOP
# ================================

SUCCESSFUL_DOWNLOADS=0
FAILED_DOWNLOADS=0

logmessage ""
logmessage "Starting RB tumor sample downloads..."
for accession in "${RB_ACCESSIONS[@]}"; do
    if download_and_convert_sample "$accession" "RB"; then
        ((SUCCESSFUL_DOWNLOADS++))
    else
        ((FAILED_DOWNLOADS++))
    fi
done

logmessage ""
logmessage "Starting control sample downloads..."
for accession in "${CONTROL_ACCESSIONS[@]}"; do
    if download_and_convert_sample "$accession" "Control"; then
        ((SUCCESSFUL_DOWNLOADS++))
    else
        ((FAILED_DOWNLOADS++))
    fi
done

# ================================
# VALIDATION
# ================================

logmessage ""
logmessage "Validating downloads..."

RB_FILES=$(find "${RAW_DATA_DIR}/RB" -name '*.fastq.gz' 2>/dev/null | wc -l)
CONTROL_FILES=$(find "${RAW_DATA_DIR}/Control" -name '*.fastq.gz' 2>/dev/null | wc -l)

logmessage "Downloaded files:"
logmessage "  - RB samples: $RB_FILES FASTQ files"
logmessage "  - Control samples: $CONTROL_FILES FASTQ files"

# Check for empty files
EMPTY_FILES=0
for file in "${RAW_DATA_DIR}"/**/*.fastq.gz; do
    if [[ -f "$file" ]] && [[ ! -s "$file" ]]; then
        logmessage "WARNING: Empty file detected: $file"
        ((EMPTY_FILES++))
    fi
done

if [[ $EMPTY_FILES -eq 0 ]]; then
    logmessage "✓ All files have content"
else
    logmessage "✗ WARNING: Found $EMPTY_FILES empty files"
fi

# ================================
# SUMMARY
# ================================

logmessage ""
logmessage "===================================================="
logmessage "DOWNLOAD PIPELINE SUMMARY"
logmessage "===================================================="
logmessage "Successful downloads: $SUCCESSFUL_DOWNLOADS"
logmessage "Failed downloads: $FAILED_DOWNLOADS"
logmessage "Download location: $RAW_DATA_DIR"
logmessage "Log file: $DOWNLOAD_LOG"
logmessage ""
logmessage "NEXT STEPS:"
logmessage "1. Review any failed downloads in the log file"
logmessage "2. Run Script 2 for Salmon quantification"
logmessage "3. Run Script 3 for Boruta feature selection"
logmessage ""
logmessage "Pipeline completed at $(date)"
logmessage "===================================================="

exit 0
