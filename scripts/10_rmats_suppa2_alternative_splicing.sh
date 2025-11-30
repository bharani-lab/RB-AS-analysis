#!/bin/bash

# Alternative Splicing Detection with rMATS and SUPPA2
# Chapter 2: Identifies differential alternative splicing events
# rMATS: Rapid analysis of multi-exon alternative splicing
# SUPPA2: Splicing Unified Pipeline for Probabilistic Assessment v2
# Author: RB-AS-analysis
# Date: 2025-11-30

# ============================================================================
# CONFIGURATION
# ============================================================================

# Output directory for rMATS results
RMATS_OUTDIR="./rMATS_results"

# SUPPA2 configuration
SUPPA2_OUTDIR="./SUPPA2_results"

# BAM file lists
BAM_LIST_RB="./bam_files/retinoblastoma_bams.txt"
BAM_LIST_CONTROL="./bam_files/control_bams.txt"

# GTF annotation file
GTF_FILE="/home/aslam/Documents/GRCh.100/Homo_sapiens.GRCh38.100.gtf"

# Reference transcriptome for SUPPA2
TRANSCRIPTOME_GTF="/home/aslam/Documents/reference/GRCh38.transcriptome.gtf"

# Python/Perl paths
PYTHON3="python3"
RMATS_PATH="/opt/rMATS_turbo_v4_1_2/rmats.py"
SUPPA2_PATH="/opt/SUPPA/suppa.py"

# Number of threads
THREADS=16

# ============================================================================
# STEP 1: rMATS Alternative Splicing Analysis
# ============================================================================

echo "=== STEP 1: rMATS Alternative Splicing Detection ==="
mkdir -v -p ${RMATS_OUTDIR}

${PYTHON3} ${RMATS_PATH} \
    --b1 ${BAM_LIST_RB} \
    --b2 ${BAM_LIST_CONTROL} \
    --gtf ${GTF_FILE} \
    --od ${RMATS_OUTDIR} \
    --tmp ${RMATS_OUTDIR}/tmp \
    -t paired \
    --readLength 101 \
    --nthread ${THREADS} \
    --statoff

echo "rMATS analysis complete. Output: ${RMATS_OUTDIR}"

# ============================================================================
# STEP 2: Extract and Process rMATS Output
# ============================================================================

echo "=== STEP 2: Processing rMATS Results ==="

# Generate summary statistics
echo "Generating rMATS summary..."
echo "Event Type,Count" > ${RMATS_OUTDIR}/rmats_summary.csv

for event_file in ${RMATS_OUTDIR}/*.MATS.JC.txt; do
    event_type=$(basename ${event_file} .MATS.JC.txt)
    count=$(tail -n +2 ${event_file} | wc -l)
    echo "${event_type},${count}" >> ${RMATS_OUTDIR}/rmats_summary.csv
done

echo "rMATS summary generated"

# ============================================================================
# STEP 3: SUPPA2 Alternative Splicing Event Prediction
# ============================================================================

echo "=== STEP 3: SUPPA2 Alternative Splicing Detection ==="
mkdir -v -p ${SUPPA2_OUTDIR}

# Generate ioe (events) file
echo "Generating ioe annotation file..."
${PYTHON3} ${SUPPA2_PATH} generateEvents \
    -i ${TRANSCRIPTOME_GTF} \
    -o ${SUPPA2_OUTDIR}/events \
    -f ioe

echo "ioe annotation file generated"

# ============================================================================
# STEP 4: Calculate PSI (Percent Spliced In) values from abundance files
# ============================================================================

echo "=== STEP 4: PSI Value Calculation ==="

# Process Salmon/Salmon transcript abundance files
for tpm_file in ./salmon_quant/*/quant.sf; do
    sample=$(basename $(dirname ${tpm_file}))
    echo "Processing ${sample} for PSI calculation..."
    
    ${PYTHON3} ${SUPPA2_PATH} psiPerEvent \
        -i ${SUPPA2_OUTDIR}/events_SE_strict.ioe \
        -e ${tpm_file} \
        -o ${SUPPA2_OUTDIR}/${sample}_psi
done

echo "PSI calculations complete"

# ============================================================================
# STEP 5: Differential Alternative Splicing Analysis
# ============================================================================

echo "=== STEP 5: Differential Alternative Splicing (DAS) Analysis ==="

# Create PSI matrix from individual samples
echo "Creating PSI matrix..."

# Merge PSI files into single matrix
paste ${SUPPA2_OUTDIR}/*_psi.psi > ${SUPPA2_OUTDIR}/psi_matrix_merged.txt

# Calculate statistics for differential ASE
echo "Calculating DAS statistics (t-test, FDR correction)..."

# Generate sample comparison
${PYTHON3} ${SUPPA2_PATH} diffSplice \
    -i ${SUPPA2_OUTDIR}/events_SE_strict.ioe \
    -e1 ${SUPPA2_OUTDIR}/RB_psi_joined.psi \
    -e2 ${SUPPA2_OUTDIR}/Control_psi_joined.psi \
    -o ${SUPPA2_OUTDIR}/RB_vs_Control

echo "Differential splicing analysis complete"

# ============================================================================
# STEP 6: Filter Significant DAS Events
# ============================================================================

echo "=== STEP 6: Filtering Significant DAS Events ==="

# Filter by p-value and PSI change threshold
awk 'NR==1 || ($NF < 0.05 && (sqrt(($6-$7)^2) >= 0.10))' \
    ${SUPPA2_OUTDIR}/RB_vs_Control.diffSplice.csv > \
    ${SUPPA2_OUTDIR}/significant_das_events.csv

echo "Filtered $(tail -n +2 ${SUPPA2_OUTDIR}/significant_das_events.csv | wc -l) significant DAS events"

# ============================================================================
# STEP 7: Generate Summary Report
# ============================================================================

echo "=== STEP 7: Generating Summary Report ==="

cat > ${SUPPA2_OUTDIR}/ANALYSIS_SUMMARY.txt << 'EOF'
Alternative Splicing Detection Pipeline Summary
===============================================

Analysis Date: $(date)

rMATS Analysis:
  - Compared RB tumors vs Normal retinal controls
  - Events detected: See rmats_summary.csv
  - Output directory: ${RMATS_OUTDIR}

SUPPA2 Analysis:
  - Events annotation: ioe format
  - PSI values calculated for all samples
  - Output directory: ${SUPPA2_OUTDIR}

Differential Alternative Splicing (DAS):
  - Comparison: RB tumors vs Control retina
  - Statistical test: t-test with FDR correction
  - Significance threshold: p < 0.05 and |ΔPSI| ≥ 0.10
  - Significant events: See significant_das_events.csv

Output Files:
  - rmats_summary.csv: Summary of rMATS-detected events
  - psi_matrix_merged.txt: PSI values across all samples
  - significant_das_events.csv: Filtered DAS events
  - RB_vs_Control.diffSplice.csv: Full differential splicing results

Next Steps:
  1. Validate DAS events with RT-qPCR (Script 12)
  2. Perform machine learning feature selection (Boruta - Script 03)
  3. Conduct protein structure analysis (Script 14)
EOF

echo "Summary report generated: ${SUPPA2_OUTDIR}/ANALYSIS_SUMMARY.txt"

echo ""
echo "=== PIPELINE COMPLETE ==="
echo "Alternative Splicing Detection Results:"
echo "  rMATS output: ${RMATS_OUTDIR}"
echo "  SUPPA2 output: ${SUPPA2_OUTDIR}"
echo "  Significant DAS events: ${SUPPA2_OUTDIR}/significant_das_events.csv"
echo ""
echo "NEXT: Run Script 11 for pan-cancer comparison"
