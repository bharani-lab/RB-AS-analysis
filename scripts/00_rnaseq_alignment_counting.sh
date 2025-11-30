#!/bin/bash

# RNA-seq Alignment and Counting Pipeline
# This script processes paired-end RNA-seq FASTQ files through quality trimming,
# STAR alignment, and featureCounts quantification
# Author: Adapted for RB-AS-analysis
# Date: 2025-11-30

# ============================================================================
# CONFIGURATION
# ============================================================================

# Path to STAR executable (trailing slash required)
STAR_PATH="/home/aslam/Documents/metaanalysis/STAR-master/bin/Linux_x86_64_static/STAR"

# fastp executable (leave empty if in PATH)
FASTP="fastp"

# Output directory prefix for STAR alignments (leave empty for current dir)
OUTDIR=""

# STAR genome index directory
INDEX="/home/aslam/Documents/metaanalysis/index"

# Reference GTF for STAR alignment
GTF="/home/aslam/Documents/GRCh.100/Homo_sapiens.GRCh38.100.gtf"

# Path to featureCounts executable (trailing slash required)
FCOUNTS_PATH="/home/aslam/Documents/metaanalysis/subread-2.0.1-Linux-x86_64/bin/"

# Output base directory for featureCounts results
FCOUNTDIR="/home/aslam/Documents/metaanalysis/olddata/featureCount"

# Filtered GTF file for featureCounts gene-level counting
FILTERED_GTF="/home/aslam/Documents/metaanalysis/ref/filter.GTF"

# ============================================================================
# STEP 1: Quality Trimming with fastp
# ============================================================================

echo "Starting RNA-seq pipeline..."
mkdir -v trim  # Create trim output directory

for i in *_1.fastq.gz; do
    SAMPLE=$(echo $i | sed 's/_1.fastq.gz//')
    echo "Processing ${SAMPLE}_1.fastq.gz and ${SAMPLE}_2.fastq.gz"
    mkdir -p -v trim/${SAMPLE}
    
    # Run fastp with 20 threads for quality trimming
    ${FASTP} \
        -i ${SAMPLE}_1.fastq.gz \
        -I ${SAMPLE}_2.fastq.gz \
        -o trim/${SAMPLE}/${SAMPLE}_1.fastq.gz \
        -O trim/${SAMPLE}/${SAMPLE}_2.fastq.gz \
        -R "${SAMPLE}" \
        -w 20
done

# ============================================================================
# STEP 2: STAR Alignment and featureCounts Quantification
# ============================================================================

cd trim

for i in *_1.fastq.gz; do
    SAMPLE=$(echo $i | sed 's/_1.fastq.gz//')
    echo "Aligning ${SAMPLE}_1.fastq.gz and ${SAMPLE}_2.fastq.gz"
    mkdir -v ${OUTDIR}${SAMPLE}
    mkdir -v ${FCOUNTDIR}${SAMPLE}
    
    # STAR 2-pass alignment with transcriptome/gene quantification
    ${STAR_PATH} \
        --runThreadN 16 \
        --genomeDir ${INDEX} \
        --sjdbGTFfile ${GTF} \
        --readFilesIn ${SAMPLE}_1.fastq.gz ${SAMPLE}_2.fastq.gz \
        --sjdbOverhang 101 \
        --readFilesCommand zcat \
        --outSAMtype BAM SortedByCoordinate \
        --outFileNamePrefix ${OUTDIR}${SAMPLE}/${SAMPLE}_ \
        --quantMode TranscriptomeSAM GeneCounts \
        --outReadsUnmapped Fastx \
        --twopassMode Basic \
        --outFilterMultimapNmax 1 \
        --outSAMstrandField intronMotif
    
    # featureCounts: gene-level count matrix from BAM files
    ${FCOUNTS_PATH}featureCounts \
        -T 30 \
        -p \
        -t exon \
        -g gene_id \
        -a ${FILTERED_GTF} \
        -o ${FCOUNTDIR}${SAMPLE}/${SAMPLE}.csv \
        ${OUTDIR}${SAMPLE}/${SAMPLE}_Aligned.sortedByCoord.out.bam
done

# ============================================================================
# STEP 3: Merge Count Matrices
# ============================================================================

cd ${FCOUNTDIR}
mkdir -v -p combine
cp *.csv combine/

cd combine

# Extract count column (column 7) from each CSV, skip header
ls -1 *.csv | parallel 'cat {} | sed 1d | cut -f7' > .clean.txt

# Create gene ID list from first CSV (column 1)
ls -1 *.csv | head -1 | xargs cut -f1 > genes.txt

# Merge into single count matrix: genes + all sample counts
paste genes.txt .clean.txt > output.txt

# Optional: Full CSV merge (all columns preserved)
paste *.csv > combine.csv

echo "Pipeline complete! Check output.txt for final count matrix."
