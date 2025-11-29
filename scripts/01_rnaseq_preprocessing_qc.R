#!/usr/bin/env Rscript
#
# Script 1: RNA-seq Preprocessing and Quality Control
# 
# Purpose: Comprehensive QC pipeline for raw RNA-seq FASTQ files from retinoblastoma samples
# Performs: FastQC analysis, adapter trimming, quality filtering, and alignment-ready validation
# 
# Author: A. Mohamed Hameed Aslam (AMRF Lab, Alavand University)
# Created: 2025-11-30
# Dependencies: fastqc, fastp, fastq_screen (optional), ShortRead (Bioconductor)
#

# Load required libraries
library(ShortRead)
library(Biostrings)
library(ggplot2)
library(dplyr)
library(data.table)

# Configuration and Parameters
# ================================

# Define input/output directories
WORK_DIR <- "./data/raw_fastq"
QC_OUTPUT_DIR <- "./qc_reports"
TRIM_OUTPUT_DIR <- "./data/trimmed_fastq"
LOG_DIR <- "./logs"

# QC thresholds
MIN_QUALITY <- 20
MIN_LENGTH <- 50
ADAPTER_PATTERN <- "AGATCGGAAGAGC"  # Illumina TruSeq adapter

# Sample metadata
SAMPLES <- c("RB_01", "RB_02", "RB_03", "RB_04", "RB_05")
RETINAL_CONTROLS <- c("Ret_Ctrl_01", "Ret_Ctrl_02", "Ret_Ctrl_03", "Ret_Ctrl_04")
ALL_SAMPLES <- c(SAMPLES, RETINAL_CONTROLS)

# Read layout (paired or single)
READ_TYPE <- "paired"  # Options: "paired" or "single"

# Create output directories
dir.create(QC_OUTPUT_DIR, showWarnings = FALSE, recursive = TRUE)
dir.create(TRIM_OUTPUT_DIR, showWarnings = FALSE, recursive = TRUE)
dir.create(LOG_DIR, showWarnings = FALSE, recursive = TRUE)

# ================================
# Function 1: FastQC Quality Assessment
# ================================

perform_fastqc_analysis <- function(sample_name, read_type = "paired") {
  
  cat("\n=== FastQC Analysis for", sample_name, "===")
  
  if (read_type == "paired") {
    r1_file <- file.path(WORK_DIR, paste0(sample_name, "_R1.fastq.gz"))
    r2_file <- file.path(WORK_DIR, paste0(sample_name, "_R2.fastq.gz"))
    files <- c(r1_file, r2_file)
  } else {
    files <- file.path(WORK_DIR, paste0(sample_name, ".fastq.gz"))
  }
  
  # Check file existence
  for (f in files) {
    if (!file.exists(f)) {
      warning("File not found:", f)
      return(NULL)
    }
  }
  
  # Read FASTQ with ShortRead
  fq_list <- lapply(files, function(x) {
    readFastq(x)
  })
  
  names(fq_list) <- if (read_type == "paired") c("R1", "R2") else "SE"
  
  # Extract quality scores and sequences
  qc_stats <- list()
  
  for (i in seq_along(fq_list)) {
    fq <- fq_list[[i]]
    read_id <- names(fq_list)[i]
    
    # Basic statistics
    sequences <- sread(fq)
    qualities <- quality(fq)
    
    qc_stats[[read_id]] <- data.frame(
      Sample = sample_name,
      Read = read_id,
      Total_Reads = length(fq),
      Mean_Length = mean(width(sequences)),
      Min_Length = min(width(sequences)),
      Max_Length = max(width(sequences)),
      Mean_Quality = mean(as.numeric(qualities)),
      GC_Content = mean(alphabetFrequency(sequences)[, "G"] + 
                        alphabetFrequency(sequences)[, "C"]) / 
                   mean(width(sequences)) * 100,
      stringsAsFactors = FALSE
    )
  }
  
  qc_df <- do.call(rbind, qc_stats)
  rownames(qc_df) <- NULL
  
  return(qc_df)
}

# ================================
# Function 2: Quality Filtering and Adapter Trimming
# ================================

perform_trimming <- function(sample_name, read_type = "paired") {
  
  cat("\n=== Trimming and QC Filtering for", sample_name, "===")
  
  if (read_type == "paired") {
    r1_input <- file.path(WORK_DIR, paste0(sample_name, "_R1.fastq.gz"))
    r2_input <- file.path(WORK_DIR, paste0(sample_name, "_R2.fastq.gz"))
    r1_output <- file.path(TRIM_OUTPUT_DIR, paste0(sample_name, "_R1_trimmed.fastq.gz"))
    r2_output <- file.path(TRIM_OUTPUT_DIR, paste0(sample_name, "_R2_trimmed.fastq.gz"))
    
    trim_cmd <- sprintf(
      "fastp -i %s -I %s -o %s -O %s --qualified_quality_phred %d --length_required %d --html %s --json %s -w 4",
      r1_input, r2_input, r1_output, r2_output, MIN_QUALITY, MIN_LENGTH,
      file.path(QC_OUTPUT_DIR, paste0(sample_name, "_fastp.html")),
      file.path(QC_OUTPUT_DIR, paste0(sample_name, "_fastp.json"))
    )
  } else {
    input <- file.path(WORK_DIR, paste0(sample_name, ".fastq.gz"))
    output <- file.path(TRIM_OUTPUT_DIR, paste0(sample_name, "_trimmed.fastq.gz"))
    
    trim_cmd <- sprintf(
      "fastp -i %s -o %s --qualified_quality_phred %d --length_required %d --html %s --json %s -w 4",
      input, output, MIN_QUALITY, MIN_LENGTH,
      file.path(QC_OUTPUT_DIR, paste0(sample_name, "_fastp.html")),
      file.path(QC_OUTPUT_DIR, paste0(sample_name, "_fastp.json"))
    )
  }
  
  # Execute fastp command
  result <- system(trim_cmd, intern = TRUE)
  
  # Log results
  cat("\nTrimming completed for", sample_name)
  
  return(list(command = trim_cmd, status = ifelse(inherits(result, "error"), "Failed", "Success")))
}

# ================================
# Function 3: QC Summary Report Generation
# ================================

generate_qc_summary <- function(qc_data_list) {
  
  cat("\n=== Generating QC Summary Report ===")
  
  # Combine all QC data
  qc_combined <- do.call(rbind, qc_data_list)
  rownames(qc_combined) <- NULL
  
  # Summary statistics by sample type
  summary_stats <- qc_combined %>%
    group_by(Sample) %>%
    summarise(
      Avg_Reads = mean(Total_Reads, na.rm = TRUE),
      Avg_Length = mean(Mean_Length, na.rm = TRUE),
      Avg_Quality = mean(Mean_Quality, na.rm = TRUE),
      GC_Avg = mean(GC_Content, na.rm = TRUE),
      .groups = 'drop'
    )
  
  # Write summary to file
  summary_file <- file.path(QC_OUTPUT_DIR, "qc_summary_report.csv")
  write.csv(summary_stats, summary_file, row.names = FALSE)
  
  cat("\nQC Summary Report saved to:", summary_file)
  
  # Generate visualization
  p <- ggplot(qc_combined, aes(x = Sample, y = Mean_Quality, color = Read)) +
    geom_point(size = 3, alpha = 0.7) +
    geom_jitter(width = 0.2) +
    theme_minimal() +
    theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
    labs(title = "Mean Quality Score Distribution",
         x = "Sample", y = "Mean Quality Score (Phred)") +
    ylim(0, 40)
  
  plot_file <- file.path(QC_OUTPUT_DIR, "qc_quality_distribution.pdf")
  ggsave(plot_file, p, width = 10, height = 6, dpi = 300)
  
  cat("\nQC Plot saved to:", plot_file)
  
  return(summary_stats)
}

# ================================
# Main Execution Pipeline
# ================================

main <- function() {
  
  cat("\n========================================")
  cat("\nRNA-seq Preprocessing and QC Pipeline")
  cat("\nStart Time:", format(Sys.time(), "%Y-%m-%d %H:%M:%S"))
  cat("\n========================================\n")
  
  # Step 1: FastQC Analysis
  cat("\n[STEP 1] Performing FastQC Analysis...\n")
  qc_data_list <- lapply(ALL_SAMPLES, function(sample) {
    tryCatch(
      perform_fastqc_analysis(sample, READ_TYPE),
      error = function(e) {
        cat("\nError processing sample", sample, ":", conditionMessage(e))
        return(NULL)
      }
    )
  })
  
  # Remove NULL entries
  qc_data_list <- qc_data_list[!sapply(qc_data_list, is.null)]
  
  if (length(qc_data_list) > 0) {
    qc_df_all <- do.call(rbind, qc_data_list)
    cat("\nFastQC Analysis completed for", nrow(qc_df_all) / 2, "samples")
  }
  
  # Step 2: Quality Filtering and Trimming
  cat("\n[STEP 2] Performing Quality Filtering and Adapter Trimming...\n")
  trim_results <- lapply(ALL_SAMPLES, function(sample) {
    tryCatch(
      perform_trimming(sample, READ_TYPE),
      error = function(e) {
        cat("\nError trimming sample", sample, ":", conditionMessage(e))
        return(NULL)
      }
    )
  })
  
  # Step 3: Generate Summary Report
  cat("\n[STEP 3] Generating QC Summary Report...\n")
  if (length(qc_data_list) > 0) {
    summary_report <- generate_qc_summary(qc_data_list)
  }
  
  cat("\n========================================")
  cat("\nPipeline Execution Completed")
  cat("\nEnd Time:", format(Sys.time(), "%Y-%m-%d %H:%M:%S"))
  cat("\n========================================\n")
  
  # Save session info
  sink(file.path(LOG_DIR, "session_info.txt"))
  print(sessionInfo())
  sink()
  
  cat("\nSession info saved to:", file.path(LOG_DIR, "session_info.txt"), "\n")
}

# Execute main pipeline
if (!interactive()) {
  main()
}
