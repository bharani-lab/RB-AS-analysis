#!/usr/bin/env Rscript
# Script 12: RT-qPCR Validation Analysis
# Chapter 3: Validates computationally predicted DAS events
# Analyzes RT-qPCR Ct values for RB vs control samples

require(ggplot2)
require(dplyr)
require(tidyr)

# Load RT-qPCR data
rtqpcr_data <- read.csv("./rtqpcr_results/ct_values.csv", row.names = 1)

# Calculate relative expression (2-ΔCt method)
relative_expr <- data.frame(
  Gene = rownames(rtqpcr_data),
  RB_mean = rowMeans(rtqpcr_data[,grep("RB", colnames(rtqpcr_data))]),
  Control_mean = rowMeans(rtqpcr_data[,grep("Control", colnames(rtqpcr_data))])
)

relative_expr$log2FC <- relative_expr$RB_mean - relative_expr$Control_mean
relative_expr$fold_change <- 2^(-relative_expr$log2FC)

# Perform t-tests
for(i in 1:nrow(rtqpcr_data)) {
  rb_vals <- as.numeric(rtqpcr_data[i, grep("RB", colnames(rtqpcr_data))])
  ctrl_vals <- as.numeric(rtqpcr_data[i, grep("Control", colnames(rtqpcr_data))])
  test <- t.test(rb_vals, ctrl_vals)
  relative_expr$p_value[i] <- test$p.value
}

relative_expr$FDR <- p.adjust(relative_expr$p_value, method = "BH")

# Visualization
png("./rtqpcr_results/validation_plot.png", width = 1200, height = 600)
p <- ggplot(relative_expr, aes(x = Gene, y = log2FC)) +
  geom_bar(stat = "identity", fill = "steelblue") +
  geom_hline(yintercept = 0, linetype = "dashed") +
  theme_minimal() +
  labs(title = "RT-qPCR Validation: Log2 Fold Change (RB vs Control)",
       y = "Log2 Fold Change", x = "Gene")
print(p)
dev.off()

# Save results
write.csv(relative_expr, "./rtqpcr_results/validation_results.csv", row.names = FALSE)
cat("\nRT-qPCR Validation Complete!\n")
print(relative_expr)
