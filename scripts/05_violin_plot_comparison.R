#!/usr/bin/env Rscript

# Script 5: Violin Plot Comparison Across Sample Groups
# Purpose: Generate violin plots for RB alternative splicing events across sample conditions
# Author: Bharani Lab
# Date: 2024

# Load required libraries
require(ggplot2)
require(ggpubr)
require(dplyr)
require(tidyr)
require(data.table)
require(scales)

# Configuration
INPUT_FILE <- "data/confirmed_rb_features.RDS"
SAMPLE_METADATA <- "data/sample_metadata.csv"
OUTPUT_DIR <- "results/05_violin_plots"

# Create output directory if it doesn't exist
if (!dir.exists(OUTPUT_DIR)) {
  dir.create(OUTPUT_DIR, recursive = TRUE)
}

# ============================================================================
# 1. LOAD DATA
# ============================================================================
logmessage("Loading RB feature data...")
confirmed_features <- readRDS(INPUT_FILE)

logmessage("Loading sample metadata...")
sample_meta <- fread(SAMPLE_METADATA)

logmessage(sprintf("Loaded %d features across %d samples", nrow(confirmed_features), ncol(confirmed_features)))

# ============================================================================
# 2. PREPARE DATA FOR PLOTTING
# ============================================================================
logmessage("Preparing data for visualization...")

# Extract feature matrix
feature_matrix <- confirmed_features[, -c(1:3)]
feature_names <- confirmed_features$feature_name

# Create long-format data for ggplot
long_data <- as.data.frame(feature_matrix) %>%
  tibble::rownames_to_column("Feature") %>%
  pivot_longer(cols = -Feature, names_to = "Sample", values_to = "Value")

# Merge with metadata
long_data <- long_data %>%
  left_join(sample_meta %>% select(Sample, Condition, Group), by = "Sample")

logmessage(sprintf("Data prepared: %d observations", nrow(long_data)))

# ============================================================================
# 3. GENERATE VIOLIN PLOTS - BY CONDITION
# ============================================================================
logmessage("Generating violin plots by condition...")

# Overall violin plot across conditions
p_condition <- ggplot(long_data, aes(x = Condition, y = Value, fill = Condition)) +
  geom_violin(alpha = 0.7, scale = "width") +
  geom_jitter(width = 0.2, height = 0, alpha = 0.3, size = 2) +
  stat_summary(fun = mean, geom = "point", shape = 21, size = 3, fill = "red") +
  theme_pubr() +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1),
    legend.position = "right",
    plot.title = element_text(hjust = 0.5, size = 14, face = "bold")
  ) +
  labs(
    title = "Distribution of RB-Specific Features by Sample Condition",
    x = "Condition",
    y = "Feature Value (Normalized)",
    fill = "Condition"
  ) +
  stat_compare_means(method = "kruskal.test", label.y = max(long_data$Value, na.rm = TRUE))

ggsave(file.path(OUTPUT_DIR, "violin_plot_condition.pdf"), p_condition, width = 10, height = 7)
logmessage("Condition violin plot saved")

# ============================================================================
# 4. GENERATE VIOLIN PLOTS - BY GROUP
# ============================================================================
logmessage("Generating violin plots by group...")

p_group <- ggplot(long_data, aes(x = Group, y = Value, fill = Group)) +
  geom_violin(alpha = 0.7, scale = "width") +
  geom_jitter(width = 0.2, height = 0, alpha = 0.3, size = 2) +
  stat_summary(fun = mean, geom = "point", shape = 21, size = 3, fill = "red") +
  facet_wrap(~Condition, scales = "free_x") +
  theme_pubr() +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1),
    legend.position = "bottom",
    plot.title = element_text(hjust = 0.5, size = 14, face = "bold"),
    strip.text = element_text(size = 11, face = "bold")
  ) +
  labs(
    title = "RB Feature Distribution by Group and Condition",
    x = "Group",
    y = "Feature Value (Normalized)",
    fill = "Group"
  )

ggsave(file.path(OUTPUT_DIR, "violin_plot_group_faceted.pdf"), p_group, width = 14, height = 8)
logmessage("Group violin plot saved")

# ============================================================================
# 5. STATISTICAL COMPARISONS
# ============================================================================
logmessage("Performing statistical tests...")

# Kruskal-Wallis test for overall difference across conditions
kw_test <- kruskal.test(Value ~ Condition, data = long_data)

logmessage(sprintf("Kruskal-Wallis test: H = %.4f, p-value = %.4e", kw_test$statistic, kw_test$p.value))

# Pairwise comparisons
pairwise_comparisons <- long_data %>%
  group_by(Feature) %>%
  summarise(
    pairwise_p = list(pairwise.wilcox.test(Value, Condition, p.adjust.method = "BH")$p.value),
    .groups = "drop"
  )

# ============================================================================
# 6. SAVE STATISTICS
# ============================================================================
logmessage("Saving statistical results...")

stats_summary <- data.table(
  Test = "Kruskal-Wallis (Overall)",
  Statistic = kw_test$statistic,
  P_Value = kw_test$p.value,
  Significant = kw_test$p.value < 0.05
)

fwrite(stats_summary, file.path(OUTPUT_DIR, "statistical_tests.csv"))

# ============================================================================
# OUTPUT SUMMARY
# ============================================================================
logmessage("")
logmessage("="*70)
logmessage("VIOLIN PLOT ANALYSIS COMPLETE")
logmessage("="*70)
logmessage(sprintf("Features analyzed: %d", length(unique(long_data$Feature))))
logmessage(sprintf("Samples analyzed: %d", length(unique(long_data$Sample))))
logmessage(sprintf("Conditions: %s", paste(unique(long_data$Condition), collapse = ", ")))
logmessage("")
logmessage("Output files:")
logmessage(sprintf(" - Condition violin plot: %s", file.path(OUTPUT_DIR, "violin_plot_condition.pdf")))
logmessage(sprintf(" - Group violin plots: %s", file.path(OUTPUT_DIR, "violin_plot_group_faceted.pdf")))
logmessage(sprintf(" - Statistics: %s", file.path(OUTPUT_DIR, "statistical_tests.csv")))
logmessage("="*70)
