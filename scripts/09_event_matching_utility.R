#!/usr/bin/env Rscript

# Script 9: Alternative Splicing Event Matching Utility
# Purpose: Match and cross-reference RB alternative splicing events across datasets and databases
# Author: Bharani Lab
# Date: 2024

# Load required libraries
require(data.table)
require(dplyr)
require(stringr)
require(igraph)

# Configuration
INPUT_EVENTS <- "data/merged/merged_rb_events.RDS"
REFERENCE_DB <- "data/external/refseq_events.csv"
OUTPUT_DIR <- "data/matched_events"

if (!dir.exists(OUTPUT_DIR)) {
  dir.create(OUTPUT_DIR, recursive = TRUE)
}

# ============================================================================
# 1. LOAD EVENT DATASETS
# ============================================================================
logmessage("Loading RB event dataset...")
rb_events <- readRDS(INPUT_EVENTS)

logmessage(sprintf("Loaded %d RB events", nrow(rb_events)))

logmessage("Loading reference database...")
ref_db <- fread(REFERENCE_DB)

logmessage(sprintf("Loaded %d reference events", nrow(ref_db)))

# ============================================================================
# 2. EVENT COORDINATE STANDARDIZATION
# ============================================================================
logmessage("Standardizing event coordinates...")

# Extract coordinates from event identifiers
rb_events <- rb_events %>%
  separate_wider_delim(
    cols = event_id,
    delim = ":",
    names = c("chr", "start", "end", "strand", "type")
  ) %>%
  mutate(
    start = as.integer(start),
    end = as.integer(end)
  )

logmessage("Event coordinates extracted and standardized")

# ============================================================================
# 3. GENE-LEVEL MATCHING
# ============================================================================
logmessage("Performing gene-level event matching...")

# Merge on gene ID
matched_genes <- rb_events %>%
  inner_join(
    ref_db %>% select(gene_id, ref_event_id, ref_annotation),
    by = "gene_id"
  )

logmessage(sprintf("Gene-level matches: %d events", nrow(matched_genes)))

# ============================================================================
# 4. COORDINATE-LEVEL MATCHING
# ============================================================================
logmessage("Performing coordinate-level event matching...")

# Define tolerance for coordinate matching (bp)
TOLERANCE <- 5

coord_matches <- rb_events %>%
  select(rb_event_id = event_id, gene_id, chr, start, end, strand, type) %>%
  as.data.table() %>%
  .[, key := paste(chr, strand, type, sep = ":")]

ref_db_processed <- ref_db %>%
  select(ref_event_id, gene_id, ref_chr = chr, ref_start = start, ref_end = end, ref_strand = strand, ref_type = type) %>%
  as.data.table() %>%
  .[, key := paste(ref_chr, ref_strand, ref_type, sep = ":")]

# Merge on key and check coordinates within tolerance
coord_matched <- merge(
  coord_matches,
  ref_db_processed,
  by = c("gene_id", "key"),
  allow.cartesian = TRUE
) %>%
  filter(
    abs(start - ref_start) <= TOLERANCE &
    abs(end - ref_end) <= TOLERANCE
  ) %>%
  mutate(
    match_type = "coordinate",
    match_score = 1 - (abs(start - ref_start) + abs(end - ref_end)) / (TOLERANCE * 2)
  )

logmessage(sprintf("Coordinate-level matches: %d events", nrow(coord_matched)))

# ============================================================================
# 5. TYPE-LEVEL CONSISTENCY CHECKING
# ============================================================================
logmessage("Checking event type consistency...")

type_consistency <- coord_matched %>%
  mutate(
    type_match = type == ref_type,
    consistency_score = ifelse(type_match, 1.0, 0.5)
  )

logmessage(sprintf("Type matches: %d (%.1f%%)", 
  sum(type_consistency$type_match),
  100 * sum(type_consistency$type_match) / nrow(type_consistency)
))

# ============================================================================
# 6. GENERATE MATCH SUMMARY
# ============================================================================
logmessage("Generating match summary...")

match_summary <- type_consistency %>%
  group_by(rb_event_id) %>%
  summarise(
    n_matches = n(),
    best_ref_id = ref_event_id[which.max(match_score)],
    best_match_score = max(match_score),
    type_consistent = any(type_match),
    .groups = "drop"
  )

# Classify matches
match_summary <- match_summary %>%
  mutate(
    match_class = case_when(
      best_match_score >= 0.95 ~ "High_Confidence",
      best_match_score >= 0.80 ~ "Medium_Confidence",
      best_match_score >= 0.60 ~ "Low_Confidence",
      TRUE ~ "Unmatched"
    )
  )

logmessage("Match classification complete")

# ============================================================================
# 7. EXPORT MATCHED EVENTS
# ============================================================================
logmessage("Exporting matched events...")

fwrite(type_consistency, file.path(OUTPUT_DIR, "detailed_matches.csv"))
fwrite(match_summary, file.path(OUTPUT_DIR, "match_summary.csv"))

# ============================================================================
# 8. MATCH STATISTICS
# ============================================================================
logmessage("Calculating match statistics...")

match_stats <- data.table(
  Category = c("Total RB Events", "Matched Events", "High Confidence", "Medium Confidence", "Low Confidence", "Unmatched"),
  Count = c(
    nrow(rb_events),
    sum(match_summary$n_matches > 0),
    sum(match_summary$match_class == "High_Confidence"),
    sum(match_summary$match_class == "Medium_Confidence"),
    sum(match_summary$match_class == "Low_Confidence"),
    sum(match_summary$match_class == "Unmatched")
  )
)

fwrite(match_stats, file.path(OUTPUT_DIR, "match_statistics.csv"))

# ============================================================================
# OUTPUT SUMMARY
# ============================================================================
logmessage("")
logmessage("="*70)
logmessage("EVENT MATCHING ANALYSIS COMPLETE")
logmessage("="*70)
logmessage(sprintf("RB events analyzed: %d", nrow(rb_events)))
logmessage(sprintf("Reference events: %d", nrow(ref_db)))
logmessage(sprintf("Total matches found: %d (%.1f%%)", 
  sum(match_summary$n_matches > 0),
  100 * sum(match_summary$n_matches > 0) / nrow(rb_events)
))
logmessage(sprintf("High confidence matches: %d", sum(match_summary$match_class == "High_Confidence")))
logmessage(sprintf("Medium confidence matches: %d", sum(match_summary$match_class == "Medium_Confidence")))
logmessage(sprintf("Low confidence matches: %d", sum(match_summary$match_class == "Low_Confidence")))
logmessage("")
logmessage("Output files generated:")
logmessage(sprintf(" - Detailed matches: %s", file.path(OUTPUT_DIR, "detailed_matches.csv")))
logmessage(sprintf(" - Match summary: %s", file.path(OUTPUT_DIR, "match_summary.csv")))
logmessage(sprintf(" - Statistics: %s", file.path(OUTPUT_DIR, "match_statistics.csv")))
logmessage("="*70)
