#!/usr/bin/env Rscript

# Script 6: Box Plot Pan-Cancer Comparison
# Purpose: Generate box plots comparing RB alternative splicing features across pan-cancer datasets
# Author: Bharani Lab
# Date: 2024

# Load required libraries
require(ggplot2)
require(ggpubr)
require(dplyr)
require(tidyr)
require(data.table)
require(beeswarm)

# Configuration
INPUT_FILE <- "data/confirmed_rb_features.RDS"
CANCER_METADATA <- "data/cancer_type_mapping.csv"
OUTPUT_DIR <- "results/06_box_plots_pan_cancer"

if (!dir.exists(OUTPUT_DIR)) {
  dir.create(OUTPUT_DIR, recursive = TRUE)
}

# ============================================================================
# 1. LOAD DATA
# ============================================================================
logmessage("Loading RB feature data and cancer metadata...")
confirmed_features <- readRDS(INPUT_FILE)
cancer_meta <- fread(CANCER_METADATA)

logmessage(sprintf("Loaded %d features", nrow(confirmed_features)))

# ============================================================================
# 2. PREPARE DATA FOR BOX PLOTS
# ============================================================================
logmessage("Preparing data for box plot visualization...")

feature_matrix <- confirmed_features[, -c(1:3)]
feature_names <- confirmed_features$feature_name

# Create long-format data
long_data <- as.data.frame(feature_matrix) %>%
  tibble::rownames_to_column("Feature") %>%
  pivot_longer(cols = -Feature, names_to = "Sample", values_to = "Value")

# Merge with cancer metadata
long_data <- long_data %>%
  left_join(cancer_meta %>% select(Sample, CancerType, Tissue), by = "Sample")

logmessage(sprintf("Pan-cancer data prepared: %d observations", nrow(long_data)))
logmessage(sprintf("Cancer types included: %s", paste(unique(long_data$CancerType), collapse = ", ")))

# ============================================================================
# 3. GENERATE BOX PLOTS - BY CANCER TYPE
# ============================================================================
logmessage("Generating box plots by cancer type...")

p_cancer <- ggplot(long_data, aes(x = CancerType, y = Value, fill = CancerType)) +
  geom_boxplot(alpha = 0.6, outlier.shape = NA) +
  geom_beeswarm(aes(color = CancerType), alpha = 0.4, size = 2, cex = 3) +
  theme_pubr() +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1, size = 11),
    legend.position = "right",
    plot.title = element_text(hjust = 0.5, size = 14, face = "bold")
  ) +
  labs(
    title = "RB-Specific Features across Pan-Cancer Datasets",
    x = "Cancer Type",
    y = "Feature Value (Normalized)",
    fill = "Cancer Type",
    color = "Cancer Type"
  ) +
  stat_compare_means(method = "kruskal.test", label.y = max(long_data$Value, na.rm = TRUE))

ggsave(file.path(OUTPUT_DIR, "box_plot_pan_cancer.pdf"), p_cancer, width = 12, height = 8)
logmessage("Pan-cancer box plot saved")

# ============================================================================
# 4. GENERATE BOX PLOTS - BY TISSUE TYPE
# ============================================================================
logmessage("Generating box plots by tissue type...")

p_tissue <- ggplot(long_data, aes(x = Tissue, y = Value, fill = CancerType)) +
  geom_boxplot(alpha = 0.6, outlier.shape = NA) +
  geom_jitter(width = 0.2, height = 0, alpha = 0.3, size = 2) +
  facet_wrap(~CancerType, scales = "free_x") +
  theme_pubr() +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1),
    legend.position = "bottom",
    plot.title = element_text(hjust = 0.5, size = 14, face = "bold"),
    strip.text = element_text(size = 10, face = "bold")
  ) +
  labs(
    title = "RB Features by Tissue Type within Cancer Categories",
    x = "Tissue Type",
    y = "Feature Value (Normalized)",
    fill = "Cancer Type"
  )

ggsave(file.path(OUTPUT_DIR, "box_plot_tissue_faceted.pdf"), p_tissue, width = 14, height = 10)
logmessage("Tissue-level box plot saved")

# ============================================================================
# 5. STATISTICAL COMPARISONS
# ============================================================================
logmessage("Performing Kruskal-Wallis test...")

kw_test <- kruskal.test(Value ~ CancerType, data = long_data)

logmessage(sprintf("Kruskal-Wallis test: H = %.4f, p-value = %.4e", kw_test$statistic, kw_test$p.value))

# ============================================================================
# 6. CANCER-TYPE SPECIFIC STATISTICS
# ============================================================================
logmessage("Calculating cancer-type specific statistics...")

cancer_stats <- long_data %>%
  group_by(CancerType) %>%
  summarise(
    n_samples = n_distinct(Sample),
    mean_value = mean(Value, na.rm = TRUE),
    median_value = median(Value, na.rm = TRUE),
    sd_value = sd(Value, na.rm = TRUE),
    .groups = "drop"
  )

fwrite(as.data.table(cancer_stats), file.path(OUTPUT_DIR, "cancer_type_statistics.csv"))

logmessage("Cancer type statistics:")
print(cancer_stats)

# ============================================================================
# 7. SUMMARY STATISTICS TABLE
# ============================================================================
logmessage("Creating summary table...")

summary_table <- data.table(
  Analysis = "Pan-Cancer Box Plots",
  Test = "Kruskal-Wallis",
  Statistic = kw_test$statistic,
  P_Value = kw_test$p.value,
  Significant = kw_test$p.value < 0.05,
  N_CancerTypes = length(unique(long_data$CancerType))
)

fwrite(summary_table, file.path(OUTPUT_DIR, "analysis_summary.csv"))

# ============================================================================
# OUTPUT SUMMARY
# ============================================================================
logmessage("")
logmessage("="*70)
logmessage("PAN-CANCER BOX PLOT ANALYSIS COMPLETE")
logmessage("="*70)
logmessage(sprintf("Total features: %d", length(unique(long_data$Feature))))
logmessage(sprintf("Total samples: %d", length(unique(long_data$Sample))))
logmessage(sprintf("Cancer types analyzed: %d", length(unique(long_data$CancerType))))
logmessage(sprintf("Tissue types: %d", length(unique(long_data$Tissue))))
logmessage("")
logmessage("Output files:")
logmessage(sprintf(" - Pan-cancer plot: %s", file.path(OUTPUT_DIR, "box_plot_pan_cancer.pdf")))
logmessage(sprintf(" - Tissue faceted plot: %s", file.path(OUTPUT_DIR, "box_plot_tissue_faceted.pdf")))
logmessage(sprintf(" - Cancer statistics: %s", file.path(OUTPUT_DIR, "cancer_type_statistics.csv")))
logmessage(sprintf(" - Summary: %s", file.path(OUTPUT_DIR, "analysis_summary.csv")))
logmessage("="*70)
