#!/bin/bash
#
# SCRIPT 2: SALMON TRANSCRIPT QUANTIFICATION AND SUPPA2 PSI CALCULATION
#
# Purpose: Quantify transcript abundances using Salmon and calculate percent spliced-in (PSI) values
# Output: TPM matrices and PSI files for alternative splicing events
# Dependencies: salmon, python3, SUPPA2
#
# Author: A. Mohamed Hameed Aslam (AMRF Lab, Allagappa University)
# Date: 2025-11-30
#

set -euo pipefail

# ================================
# LOGGING FUNCTION
# ================================

logmsg() {
    local timestamp=$(date '+%Y-%m-%d %H:%M:%S')
    echo "[$timestamp] $1"
    echo "[$timestamp] $1" >> "$LOG_FILE"
}

# ================================
# CONFIGURATION
# ================================

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
PROJECT_ROOT="${SCRIPT_DIR%/scripts}"
FASTQ_INPUT_DIR="${PROJECT_ROOT}/data/trimmed_fastq"  # Output from Script 1
SALMON_OUTPUT_DIR="${PROJECT_ROOT}/results/salmon_quantification"
SUPPA_OUTPUT_DIR="${PROJECT_ROOT}/results/suppa_psi"
LOG_DIR="${PROJECT_ROOT}/logs"
LOG_FILE="${LOG_DIR}/salmon_quantification_$(date +%Y%m%d_%H%M%S).log"

mkdir -p "$SALMON_OUTPUT_DIR"
mkdir -p "$SUPPA_OUTPUT_DIR"
mkdir -p "$LOG_DIR"

# ================================
# USER-CONFIGURED PATHS
# ================================
# UPDATE THESE PATHS TO YOUR INSTALLATION:
SALMON_INDEX_DIR="/path/to/salmon/index"  # Pre-built Salmon index for your reference
SUPPA_ANNOTATION="/path/to/suppa/annotation.ioe"  # SUPPA2 isoform event file
SUPPA_SCRIPT_DIR="/path/to/SUPPA2/scripts"  # SUPPA2 Python scripts directory
THREADS=8  # Number of CPU threads to use
BOOTSTRAP_SAMPLES=30  # Bootstrap replicates for uncertainty estimation
LIBRARY_TYPE="ISR"  # Illumina Stranded Reverse

logmsg "Configuration:"
logmsg "  - Salmon Index: $SALMON_INDEX_DIR"
logmsg "  - SUPPA Annotation: $SUPPA_ANNOTATION"
logmsg "  - SUPPA Scripts: $SUPPA_SCRIPT_DIR"
logmsg "  - Input FASTQ: $FASTQ_INPUT_DIR"
logmsg "  - Output Directory: $SALMON_OUTPUT_DIR"
logmsg "  - Threads: $THREADS"

# ================================
# DEPENDENCY CHECK
# ================================

if ! command -v salmon &> /dev/null; then
    logmsg "ERROR: Salmon not found. Install with: conda install -c bioconda salmon"
    exit 1
fi

if ! command -v python3 &> /dev/null; then
    logmsg "ERROR: Python3 not found"
    exit 1
fi

logmsg "Salmon version: $(salmon --version)"

# ================================
# SALMON QUANTIFICATION FUNCTION
# ================================

quantify_sample() {
    local r1_file="$1"
    local sample_name=$(basename "$r1_file" _R1_trimmed.fastq.gz)
    local r2_file="${r1_file%_R1_trimmed.fastq.gz}_R2_trimmed.fastq.gz"
    local sample_salmon_dir="${SALMON_OUTPUT_DIR}/${sample_name}"
    
    if [[ ! -f "$r2_file" ]]; then
        logmsg "WARNING: Paired file not found for $sample_name"
        return 1
    fi
    
    logmsg "Processing sample: $sample_name"
    
    if salmon quant \
        --index "$SALMON_INDEX_DIR" \
        --libType "$LIBRARY_TYPE" \
        --mates1 "$r1_file" \
        --mates2 "$r2_file" \
        --output "$sample_salmon_dir" \
        --threads "$THREADS" \
        --validateMappings \
        --numBootstraps "$BOOTSTRAP_SAMPLES" 2>> "$LOG_FILE"; then
        logmsg "Salmon quantification successful for $sample_name"
        
        # Extract TPM values (isoform level)
        if [[ -f "${sample_salmon_dir}/quant.sf" ]]; then
            python3 "${SUPPA_SCRIPT_DIR}/multipleFieldSelection.py" \
                --input "${sample_salmon_dir}/quant.sf" \
                --key 1 \
                --field 4 \
                --output "${sample_salmon_dir}/${sample_name}_isotpm.txt" 2>> "$LOG_FILE"
            logmsg "Isoform TPM extraction successful"
        fi
        
        # Calculate PSI values
        if [[ -f "${sample_salmon_dir}/${sample_name}_isotpm.txt" ]]; then
            python3 "${SUPPA_SCRIPT_DIR}/psiPerEvent.py" \
                --ioe-file "$SUPPA_ANNOTATION" \
                --expression-file "${sample_salmon_dir}/${sample_name}_isotpm.txt" \
                --output-file "${SUPPA_OUTPUT_DIR}/${sample_name}" 2>> "$LOG_FILE"
            logmsg "PSI calculation successful for $sample_name"
        fi
        
        return 0
    else
        logmsg "ERROR: Salmon quantification failed for $sample_name"
        return 1
    fi
}

# ================================
# MAIN PROCESSING LOOP
# ================================

logmsg ""
logmsg "===================================================="
logmsg "Starting Salmon Quantification Pipeline"
logmsg "===================================================="
logmsg ""

SUCCESS_COUNT=0
FAILED_COUNT=0

# Find all R1 FASTQ files
cd "$FASTQ_INPUT_DIR" 2>/dev/null || {
    logmsg "ERROR: FASTQ input directory not found: $FASTQ_INPUT_DIR"
    exit 1
}

for r1_file in *_R1_trimmed.fastq.gz; do
    if [[ -f "$r1_file" ]]; then
        if quantify_sample "${FASTQ_INPUT_DIR}/${r1_file}"; then
            ((SUCCESS_COUNT++))
        else
            ((FAILED_COUNT++))
        fi
    fi
done

# ================================
# CREATE MERGED PSI MATRIX
# ================================

logmsg ""
logmsg "Creating consolidated PSI matrix..."

if [[ -d "$SUPPA_OUTPUT_DIR" ]] && [[ $(find "$SUPPA_OUTPUT_DIR" -name '*.psi' | wc -l) -gt 0 ]]; then
    # Python script to merge PSI files
    python3 << 'EOF'
import os
import pandas as pd
from glob import glob

suppa_dir = os.environ['SUPPA_OUTPUT_DIR']
psi_files = glob(f'{suppa_dir}/*.psi')

if psi_files:
    psi_data = []
    for psi_file in sorted(psi_files):
        df = pd.read_csv(psi_file, sep='\t', index_col=0)
        sample_name = os.path.basename(psi_file).replace('.psi', '')
        df.columns = [sample_name]
        psi_data.append(df)
    
    merged_psi = pd.concat(psi_data, axis=1)
    merged_psi.to_csv(f'{suppa_dir}/merged_psi_matrix.tsv', sep='\t')
    print(f'PSI matrix created: {suppa_dir}/merged_psi_matrix.tsv')
EOF
fi

# ================================
# SUMMARY
# ================================

logmsg ""
logmsg "===================================================="
logmsg "QUANTIFICATION PIPELINE SUMMARY"
logmsg "===================================================="
logmsg "Successful quantifications: $SUCCESS_COUNT"
logmsg "Failed quantifications: $FAILED_COUNT"
logmsg "Salmon output: $SALMON_OUTPUT_DIR"
logmsg "SUPPA PSI output: $SUPPA_OUTPUT_DIR"
logmsg "Log file: $LOG_FILE"
logmsg ""
logmsg "NEXT STEPS:"
logmsg "1. Review PSI matrix at: $SUPPA_OUTPUT_DIR/merged_psi_matrix.tsv"
logmsg "2. Run Script 3 for Boruta feature selection"
logmsg ""
logmsg "Pipeline completed at $(date)"
logmsg "===================================================="

exit 0
