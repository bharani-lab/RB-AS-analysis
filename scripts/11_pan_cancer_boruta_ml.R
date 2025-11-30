#!/usr/bin/env Rscript
# Script 11: Pan-Cancer Boruta Machine Learning
# Chapter 2: Identifies RB-specific alternative splicing events
# Uses Boruta algorithm to distinguish RB from other 33 TCGA cancer types
# Author: RB-AS-analysis
# Date: 2025-11-30

# ============================================================================
# SETUP & LIBRARIES
# ============================================================================

require(Boruta)
require(randomForest)
require(ggplot2)
require(dplyr)

set.seed(42)

# Configuration
RB_DAS_FILE <- "./SUPPA2_results/significant_das_events.csv"
PAN_CANCER_PSI_FILE <- "./pan_cancer_data/tcga_psi_matrix.csv"
OUTPUT_DIR <- "./Boruta_ML_results"

dir.create(OUTPUT_DIR, showWarnings = FALSE)

# ============================================================================
# STEP 1: Load DAS Events & Pan-Cancer Data
# ============================================================================

cat("\n=== Loading RB DAS events and pan-cancer data ===\n")

rb_das <- read.csv(RB_DAS_FILE, row.names = 1)
pan_cancer_psi <- read.csv(PAN_CANCER_PSI_FILE, row.names = 1)

# Create binary classification: RB=1, Others=0
rb_samples <- colnames(rb_das)
other_samples <- setdiff(colnames(pan_cancer_psi), rb_samples)

cat(paste("RB samples:", length(rb_samples), "\n"))
cat(paste("Other cancer samples:", length(other_samples), "\n"))

# Merge datasets
merged_data <- pan_cancer_psi[, c(rb_samples, other_samples)]
target <- c(rep(1, length(rb_samples)), rep(0, length(other_samples)))

# ============================================================================
# STEP 2: Boruta Feature Selection
# ============================================================================

cat("\n=== Running Boruta Algorithm (1000 iterations) ===\n")

boruta_result <- Boruta(as.factor(target) ~ ., 
                        data = t(merged_data),
                        maxRuns = 1000,
                        pValue = 0.01,
                        doTrace = 2)

cat("\nBoruta Analysis Complete!\n")

# Extract confirmed features
confirmed_features <- getSelectedAttributes(boruta_result, withTentative = FALSE)
confirmed_tentative <- getSelectedAttributes(boruta_result, withTentative = TRUE)

cat(paste("\nConfirmed RB-specific DAS events:", length(confirmed_features), "\n"))
cat(paste("Confirmed + Tentative:", length(confirmed_tentative), "\n"))

# ============================================================================
# STEP 3: Generate Importance Plot
# ============================================================================

cat("\n=== Generating Importance Plots ===\n")

png(file.path(OUTPUT_DIR, "boruta_importance.png"), width = 1200, height = 800)
plot(boruta_result, las = 2, main = "Boruta Feature Importance: RB vs Other Cancers")
dev.off()

# ============================================================================
# STEP 4: Random Forest Classification Performance
# ============================================================================

cat("\n=== Evaluating Random Forest Classification ===\n")

# Train RF on confirmed features
rf_data <- t(merged_data[confirmed_features, ])
rf_model <- randomForest(as.factor(target) ~ .,
                         data = as.data.frame(rf_data),
                         ntree = 500,
                         importance = TRUE)

print(rf_model)

# ============================================================================
# STEP 5: Cross-Validation
# ============================================================================

cat("\n=== 10-Fold Cross-Validation ===\n")

folds <- cut(seq_along(target), breaks = 10, labels = FALSE)
accuracy_scores <- numeric(10)

for (fold in 1:10) {
  test_idx <- which(folds == fold)
  train_data <- rf_data[-test_idx, ]
  test_data <- rf_data[test_idx, ]
  train_target <- target[-test_idx]
  test_target <- target[test_idx]
  
  fold_rf <- randomForest(as.factor(train_target) ~ .,
                          data = as.data.frame(train_data),
                          ntree = 500)
  
  predictions <- predict(fold_rf, as.data.frame(test_data))
  accuracy_scores[fold] <- mean(predictions == test_target)
}

cat(paste("\nMean CV Accuracy:", round(mean(accuracy_scores), 4), "\n"))
cat(paste("SD:", round(sd(accuracy_scores), 4), "\n"))

# ============================================================================
# STEP 6: Output Summary & Confirmed Features
# ============================================================================

cat("\n=== Saving Results ===\n")

# Save confirmed features
confirmed_df <- data.frame(
  Event = confirmed_features,
  Importance = rf_model$importance[confirmed_features, "MeanDecreaseGini"]
) %>% arrange(desc(Importance))

write.csv(confirmed_df, 
          file.path(OUTPUT_DIR, "confirmed_rb_specific_events.csv"),
          row.names = FALSE)

# Summary statistics
summary_text <- sprintf(
  "RB-Specific Alternative Splicing ML Analysis\n%s\n\nBoruta Algorithm Results:\n- Confirmed RB-specific events: %d\n- Confirmed + Tentative: %d\n- Random Forest CV Accuracy: %.4f (SD: %.4f)\n\nNext Steps:\n1. Validate top events with RT-qPCR (Script 12)\n2. Conduct protein structure analysis (Script 14)\n3. Proteomics validation (Script 13)\n",
  format(Sys.time(), "%Y-%m-%d %H:%M"),
  length(confirmed_features),
  length(confirmed_tentative),
  mean(accuracy_scores),
  sd(accuracy_scores)
)

writeLines(summary_text, file.path(OUTPUT_DIR, "ANALYSIS_SUMMARY.txt"))

cat("\n=== PIPELINE COMPLETE ===\n")
cat(sprintf("Results saved to: %s\n", OUTPUT_DIR))
cat(sprintf("Confirmed RB-specific events: %d\n", length(confirmed_features)))
cat(sprintf("Random Forest CV Accuracy: %.4f\n\n", mean(accuracy_scores)))
