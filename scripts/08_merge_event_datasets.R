#!/usr/bin/env Rscript

# Script 8: Merge Alternative Splicing Event Datasets
# Purpose: Consolidate multiple RB alternative splicing event datasets and prepare unified data structure
# Author: Bharani Lab
# Date: 2024

# Load required libraries
require(data.table)
require(dplyr)
require(tidyr)

# Configuration
INPUT_DIR <- "data/raw_events"
OUTPUT_DIR <- "data/merged"

if (!dir.exists(OUTPUT_DIR)) {
  dir.create(OUTPUT_DIR, recursive = TRUE)
}

# ============================================================================
# 1. DISCOVER INPUT FILES
# ============================================================================
logmessage("Discovering event dataset files...")

event_files <- list.files(INPUT_DIR, pattern = "\\.csv$", full.names = TRUE)
logmessage(sprintf("Found %d event files", length(event_files)))

# ============================================================================
# 2. LOAD AND MERGE DATASETS
# ============================================================================
logmessage("Loading and merging event datasets...")

all_events <- data.table()

for (file in event_files) {
  logmessage(sprintf("Processing: %s", basename(file)))
  
  events <- fread(file)
  
  # Add source file information
  events$source_file <- basename(file)
  
  # Bind to combined dataset
  all_events <- rbindlist(list(all_events, events), fill = TRUE, use.names = TRUE)
}

logmessage(sprintf("Total events loaded: %d", nrow(all_events)))

# ============================================================================
# 3. DATA STANDARDIZATION
# ============================================================================
logmessage("Standardizing event data...")

# Ensure consistent column types
all_events$event_id <- as.character(all_events$event_id)
all_events$gene_id <- as.character(all_events$gene_id)
all_events$event_type <- as.character(all_events$event_type)

# Identify and handle duplicates
dup_count <- sum(duplicated(all_events$event_id))
logmessage(sprintf("Duplicate events found: %d", dup_count))

# Remove duplicates, keeping first occurrence
all_events <- all_events[!duplicated(all_events$event_id)]

logmessage(sprintf("Final unique events: %d", nrow(all_events)))

# ============================================================================
# 4. EVENT TYPE SUMMARY
# ============================================================================
logmessage("Summarizing event types...")

event_summary <- all_events %>%
  group_by(event_type) %>%
  summarise(
    count = n(),
    unique_genes = n_distinct(gene_id),
    .groups = "drop"
  ) %>%
  arrange(desc(count))

logmessage("Event type distribution:")
print(event_summary)

# ============================================================================
# 5. GENE COVERAGE ANALYSIS
# ============================================================================
logmessage("Analyzing gene coverage...")

gene_coverage <- all_events %>%
  group_by(gene_id) %>%
  summarise(
    total_events = n(),
    event_types = n_distinct(event_type),
    .groups = "drop"
  )

logmessage(sprintf("Total genes with events: %d", nrow(gene_coverage)))
logmessage(sprintf("Avg events per gene: %.2f", mean(gene_coverage$total_events)))

# ============================================================================
# 6. EXPORT MERGED DATASET
# ============================================================================
logmessage("Exporting merged dataset...")

fwrite(all_events, file.path(OUTPUT_DIR, "merged_rb_events.csv"))
fwrite(as.data.table(event_summary), file.path(OUTPUT_DIR, "event_type_summary.csv"))
fwrite(gene_coverage, file.path(OUTPUT_DIR, "gene_coverage.csv"))

logmessage("Datasets exported successfully")

# ============================================================================
# 7. SAVE RDS FORMAT
# ============================================================================
logmessage("Saving RDS format for R analysis...")

saveRDS(all_events, file.path(OUTPUT_DIR, "merged_rb_events.RDS"))

logmessage("RDS file saved")

# ============================================================================
# OUTPUT SUMMARY
# ============================================================================
logmessage("")
logmessage("="*70)
logmessage("EVENT DATASET MERGE COMPLETE")
logmessage("="*70)
logmessage(sprintf("Input files processed: %d", length(event_files)))
logmessage(sprintf("Total unique events: %d", nrow(all_events)))
logmessage(sprintf("Total genes covered: %d", nrow(gene_coverage)))
logmessage(sprintf("Event types: %d", nrow(event_summary)))
logmessage("")
logmessage(sprintf("Merged data saved to: %s", OUTPUT_DIR))
logmessage("="*70)
