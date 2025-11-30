#!/usr/bin/env Rscript
# Script 13: Mass Spectrometry Proteomics Analysis
# Chapter 3: Validates DAS events at protein level
require(dplyr)
require(ggplot2)

# Load MaxQuant proteomics data
peptides <- read.csv("./proteomics_data/peptides.csv", row.names = 1)

# Extract isoform-specific peptides
isoform_peps <- peptides %>%
  filter(Isoform_specific == TRUE) %>%
  group_by(Protein) %>%
  summarise(
    n_unique_peptides = n_distinct(Sequence),
    total_psms = sum(PSM_count),
    seq_coverage = mean(Coverage)
  )

cat(sprintf("\nProteomics Summary:\n"))
cat(sprintf("Proteins with detected isoforms: %d\n", nrow(isoform_peps)))
cat(sprintf("Avg sequence coverage: %.2f%%\n", mean(isoform_peps$seq_coverage)))

# Save results
summary_stats <- data.frame(
  Protein = isoform_peps$Protein,
  Unique_Peptides = isoform_peps$n_unique_peptides,
  Total_PSMs = isoform_peps$total_psms,
  Sequence_Coverage = isoform_peps$seq_coverage
)

write.csv(summary_stats, "./proteomics_results/isoform_summary.csv", row.names = FALSE)
cat("\nMass Spectrometry Analysis Complete!\n")
