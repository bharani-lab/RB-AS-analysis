#!/usr/bin/env Rscript
#
# SCRIPT 3: BORUTA MACHINE LEARNING FEATURE SELECTION
#
# Purpose: Identify RB-specific alternative splicing events using Boruta algorithm
# Input: PSI matrix from Script 2 (Salmon/SUPPA2 quantification)
# Output: Confirmed RB-specific features, importance rankings, and plots
# Dependencies: readr, tidyverse, Boruta, impute, randomForest
#
# Author: A. Mohamed Hameed Aslam (AMRF Lab, Allagappa University)
# Date: 2025-11-30
#

suppressPackageStartupMessages({
  library(readr)        # Fast data import/export
  library(tidyverse)    # Data manipulation and visualization
  library(Boruta)       # Feature selection algorithm
  library(impute)       # k-NN imputation for missing values
  library(randomForest) # Random forest backend for Boruta
})

cat("\nINFO: All required libraries loaded successfully\n")

# ================================
# CONFIGURATION
# ================================

PROJECT_DIR <- "/path/to/your/project"  # CUSTOMIZE: Replace with your project root
OUTPUT_DIR <- file.path(PROJECT_DIR, "results/feature_selection")
DATA_DIR <- file.path(PROJECT_DIR, "data/processed")
LOG_DIR <- file.path(PROJECT_DIR, "logs")

dir.create(OUTPUT_DIR, showWarnings = FALSE, recursive = TRUE)
dir.create(LOG_DIR, showWarnings = FALSE, recursive = TRUE)

LOG_FILE <- file.path(LOG_DIR, paste0("boruta_", Sys.Date(), ".log"))

# ================================
# LOGGING FUNCTION
# ================================

logmessage <- function(msg, level = "INFO") {
  timestamp <- format(Sys.time(), "%Y-%m-%d %H:%M:%S")
  logtext <- sprintf("[%s] %s: %s", timestamp, level, msg)
  cat(logtext, "\n")
  cat(logtext, "\n", file = LOG_FILE, append = TRUE)
}

logmessage("Starting Boruta Feature Selection Pipeline")
logmessage(paste("Project directory:", PROJECT_DIR))
logmessage(paste("Output directory:", OUTPUT_DIR))

# ================================
# LOAD PSI DATA
# ================================

logmessage("Loading PSI matrix...")
psi_data <- read_csv(file.path(DATA_DIR, "merged_psi_matrix.tsv"),
                      col_types = cols(.default = col_double(),
                                      EventID = col_character()))

logmessage(sprintf("Data loaded: %d events x %d samples", nrow(psi_data), ncol(psi_data) - 1))

# Extract event IDs and PSI values
event_ids <- psi_data$EventID
psi_matrix <- psi_data %>% select(-EventID) %>% as.data.frame()
rownames(psi_matrix) <- event_ids

# ================================
# SAMPLE ANNOTATIONS
# ================================

# Define sample groups (RB: 1-32, Control: 33-57)
num_rb <- 32
num_control <- 25
sample_groups <- c(rep("RB", num_rb), rep("Control", num_control))
names(sample_groups) <- colnames(psi_matrix)

logmessage(sprintf("RB samples: %d, Control samples: %d", num_rb, num_control))

# ================================
# K-NN IMPUTATION (HANDLE MISSING VALUES)
# ================================

logmessage("Performing k-NN imputation for missing values...")
psi_matrix_numeric <- as.matrix(psi_matrix)
storage.mode(psi_matrix_numeric) <- "numeric"

if (any(is.na(psi_matrix_numeric))) {
  psi_matrix_imputed <- impute.knn(psi_matrix_numeric, k = 10)$data
  logmessage("k-NN imputation completed")
} else {
  psi_matrix_imputed <- psi_matrix_numeric
  logmessage("No missing values detected")
}

# ================================
# PREPARE DATA FOR BORUTA
# ================================

logmessage("Preparing data for Boruta algorithm...")

# Create classification vector: RB = 1, Control = 0
class_vector <- as.factor(ifelse(sample_groups == "RB", "RB", "Control"))

# Create data frame for Boruta
boruta_data <- data.frame(
  Class = class_vector,
  t(psi_matrix_imputed)
)

logmessage(sprintf("Boruta input: %d samples x %d features", nrow(boruta_data), ncol(boruta_data) - 1))

# ================================
# RUN BORUTA ALGORITHM
# ================================

logmessage("Running Boruta feature selection (this may take several minutes)...")

set.seed(42)  # For reproducibility
boruta_result <- Boruta(Class ~ ., data = boruta_data, doTrace = 2, ntree = 1000)

logmessage("Boruta algorithm completed")

# ================================
# EXTRACT CONFIRMED FEATURES
# ================================

logmessage("Extracting confirmed RB-specific features...")

confirmed_features <- getSelectedAttributes(boruta_result, withTentative = FALSE)
confirmed_importance <- attStats(boruta_result)
confirmed_importance <- confirmed_importance[confirmed_features, ]
confirmed_importance <- confirmed_importance[order(confirmed_importance$medianImp, decreasing = TRUE), ]

logmessage(sprintf("Confirmed RB-specific AS events: %d", nrow(confirmed_importance)))
logmessage(sprintf("Rejected non-informative features: %d", 
                   sum(boruta_result$finalDecision == "Rejected")))

# ================================
# SAVE RESULTS
# ================================

logmessage("Saving results...")

# Save confirmed features
output_confirmed <- rownames_to_column(as.data.frame(confirmed_importance), "EventID")
write_csv(output_confirmed, file.path(OUTPUT_DIR, "rb_confirmed_features.csv"))
logmessage("Saved: rb_confirmed_features.csv")

# Save all importance statistics
output_all <- rownames_to_column(as.data.frame(attStats(boruta_result)), "EventID")
output_all <- arrange(output_all, desc(medianImp))
write_csv(output_all, file.path(OUTPUT_DIR, "importance_statistics_all.csv"))
logmessage("Saved: importance_statistics_all.csv")

# ================================
# GENERATE PLOTS
# ================================

logmessage("Creating importance plots...")

pdf(file.path(OUTPUT_DIR, "boruta_importance_plot.pdf"), width = 14, height = 10)
plot(boruta_result, cex.axis = 0.7, las = 2,
     main = "Boruta Feature Selection: RB-Specific Alternative Splicing Events",
     xlab = "Alternative Splicing Events",
     ylab = "Importance Score")
dev.off()
logmessage("Saved: boruta_importance_plot.pdf")

# ================================
# SUMMARY REPORT
# ================================

logmessage("\n" %+% strrep("=", 70))
logmessage("BORUTA FEATURE SELECTION SUMMARY REPORT")
logmessage(strrep("=", 70))
logmessage("")
logmessage("Analysis Overview:")
logmessage(sprintf("  Total features analyzed: %d", ncol(boruta_data) - 1))
logmessage(sprintf("  RB samples: %d", num_rb))
logmessage(sprintf("  Control samples: %d", num_control))
logmessage("")
logmessage("Results:")
logmessage(sprintf("  Confirmed robust biomarkers: %d (%.1f%%)",
                   nrow(confirmed_importance),
                   nrow(confirmed_importance) / (ncol(boruta_data) - 1) * 100))
logmessage(sprintf("  Rejected not informative: %d (%.1f%%)",
                   sum(boruta_result$finalDecision == "Rejected"),
                   sum(boruta_result$finalDecision == "Rejected") / (ncol(boruta_data) - 1) * 100))
logmessage(sprintf("  Tentative borderline: %d (%.1f%%)",
                   sum(boruta_result$finalDecision == "Tentative"),
                   sum(boruta_result$finalDecision == "Tentative") / (ncol(boruta_data) - 1) * 100))
logmessage("")
logmessage("Top 10 Most Important RB-Specific Features:")
for (i in 1:min(10, nrow(confirmed_importance))) {
  event_name <- rownames(confirmed_importance)[i]
  median_imp <- confirmed_importance$medianImp[i]
  logmessage(sprintf("  %2d. %s - Importance: %.4f", i, event_name, median_imp))
}
logmessage("")
logmessage("Output Files:")
logmessage(sprintf("  - %s", file.path(OUTPUT_DIR, "rb_confirmed_features.csv")))
logmessage(sprintf("  - %s", file.path(OUTPUT_DIR, "importance_statistics_all.csv")))
logmessage(sprintf("  - %s", file.path(OUTPUT_DIR, "boruta_importance_plot.pdf")))
logmessage("")
logmessage("NEXT STEPS:")
logmessage("1. Review confirmed features for RB-specific AS patterns")
logmessage("2. Perform enrichment analysis on confirmed genes")
logmessage("3. Validate top features experimentally (RT-qPCR, minigene assay)")
logmessage("4. Build classification/regression models using confirmed features")
logmessage("")
logmessage("Pipeline completed successfully!")
logmessage(strrep("=", 70))

cat("\n✓ Boruta feature selection pipeline complete!\n")
cat("  Output directory:", OUTPUT_DIR, "\n\n")
