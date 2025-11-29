#!/usr/bin/env Rscript

# Script 4: Hierarchical Clustering & Heatmap Visualization
# Purpose: Perform hierarchical clustering on RB alternative splicing data and generate publication-quality heatmaps
# Author: Bharani Lab
# Date: 2024

# Load required libraries
require(data.table)
require(ComplexHeatmap)
require(circlize)
require(cluster)
require(dendextend)

# Configuration
INPUT_FILE <- "data/confirmed_rb_features.RDS"
OUTPUT_DIR <- "results/04_clustering_heatmap"

# Create output directory if it doesn't exist
if (!dir.exists(OUTPUT_DIR)) {
  dir.create(OUTPUT_DIR, recursive = TRUE)
}

# ============================================================================
# 1. LOAD AND PREPARE DATA
# ============================================================================
logmessage("Loading confirmed RB-specific features...")
confirmed_features <- readRDS(INPUT_FILE)

logmessage(sprintf("Loaded %d features across %d samples", nrow(confirmed_features), ncol(confirmed_features)))

# Extract feature matrix (remove metadata columns)
feature_matrix <- confirmed_features[, -c(1:3)]  # Assuming first 3 cols are metadata
feature_matrix <- as.matrix(feature_matrix)

# Standardize features for clustering (z-score normalization)
feature_matrix_scaled <- scale(feature_matrix, center = TRUE, scale = TRUE)

logmessage("Data standardized for clustering")

# ============================================================================
# 2. HIERARCHICAL CLUSTERING ANALYSIS
# ============================================================================
logmessage("Performing hierarchical clustering...")

# Calculate distance matrix using Euclidean distance
dist_matrix <- dist(feature_matrix_scaled, method = "euclidean")

# Perform hierarchical clustering using Ward.D2 method
hc_result <- hclust(dist_matrix, method = "ward.D2")

logmessage(sprintf("Hierarchical clustering completed"))

# ============================================================================
# 3. DENDROGRAM GENERATION AND CUTOFF
# ============================================================================
logmessage("Generating dendrograms...")

# Cut dendrogram at optimal height (using elbow method or fixed cutoff)
optimal_clusters <- 4  # Can be optimized using silhouette analysis
cluster_assignments <- cutree(hc_result, k = optimal_clusters)

logmessage(sprintf("Data partitioned into %d clusters", optimal_clusters))

# Create dendrogram with color-coded branches
dend <- as.dendrogram(hc_result)
dend <- color_branches(dend, k = optimal_clusters, col = c("#E41A1C", "#377EB8", "#4DAF4A", "#984EA3"))

# ============================================================================
# 4. COMPLEX HEATMAP GENERATION
# ============================================================================
logmessage("Generating publication-quality heatmap...")

# Reorder data by hierarchical clustering
row_order <- hc_result$order
feature_matrix_ordered <- feature_matrix_scaled[row_order, ]

# Define color gradient for heatmap
col_fun <- colorRamp2(c(-2, 0, 2), c("#0072B2", "white", "#D55E00"))

# Create sample annotations if available
sample_anno_df <- data.frame(
  Cluster = factor(cluster_assignments[row_order])
)

sample_anno <- HeatmapAnnotation(
  Cluster = sample_anno_df$Cluster,
  col = list(Cluster = c("1" = "#E41A1C", "2" = "#377EB8", "3" = "#4DAF4A", "4" = "#984EA3")),
  annotation_label = "Cluster Assignment",
  show_annotation_name = TRUE
)

# Generate main heatmap
main_heatmap <- Heatmap(
  feature_matrix_ordered,
  name = "Scaled Expression",
  col = col_fun,
  cluster_rows = FALSE,  # Already ordered by hierarchical clustering
  cluster_columns = TRUE,
  top_annotation = sample_anno,
  row_dend_side = "left",
  column_dend_side = "top",
  show_row_names = FALSE,
  show_column_names = TRUE,
  row_title = "RB-Specific Features",
  column_title = "Samples",
  heatmap_legend_param = list(
    title = "Z-Score",
    title_position = "topcenter",
    legend_direction = "vertical"
  )
)

# Save heatmap as high-resolution PDF
pdf(file.path(OUTPUT_DIR, "hierarchical_clustering_heatmap.pdf"), width = 16, height = 12)
print(main_heatmap)
dev.off()

logmessage("Heatmap saved to PDF")

# ============================================================================
# 5. CLUSTER STATISTICS AND CHARACTERIZATION
# ============================================================================
logmessage("Calculating cluster statistics...")

# Calculate silhouette scores for cluster quality assessment
silhouette_scores <- silhouette(cluster_assignments, dist_matrix)
avg_silhouette <- mean(silhouette_scores[, 3])

logmessage(sprintf("Average silhouette score: %.4f", avg_silhouette))

# Per-cluster statistics
cluster_stats <- data.table(
  Cluster = 1:optimal_clusters,
  Size = tabulate(cluster_assignments),
  Avg_Silhouette = as.numeric(by(silhouette_scores[, 3], cluster_assignments, mean)),
  Within_Cluster_SS = as.numeric(by(dist_matrix, cluster_assignments, function(x) sum(x^2) / length(x)))
)

logmessage("Cluster Statistics:")
print(cluster_stats)

# ============================================================================
# 6. EXPORT CLUSTER ASSIGNMENTS
# ============================================================================
logmessage("Exporting cluster assignments...")

# Create assignment table
cluster_assignment_table <- data.table(
  Feature = rownames(feature_matrix)[row_order],
  Cluster = cluster_assignments[row_order],
  Silhouette_Score = silhouette_scores[row_order, 3]
)

# Save assignments
fwrite(cluster_assignment_table, file.path(OUTPUT_DIR, "cluster_assignments.csv"))

logmessage("Cluster assignments saved to CSV")

# ============================================================================
# 7. DENDROGRAM EXPORT
# ============================================================================
logmessage("Exporting dendrogram...")

pdf(file.path(OUTPUT_DIR, "dendrogram_colored.pdf"), width = 14, height = 10)
plot(dend, main = "Hierarchical Clustering Dendrogram (Ward.D2 Method)")
dev.off()

logmessage("Dendrogram saved to PDF")

# ============================================================================
# OUTPUT SUMMARY
# ============================================================================
logmessage("")
logmessage("="*70)
logmessage("HIERARCHICAL CLUSTERING ANALYSIS COMPLETE")
logmessage("="*70)
logmessage(sprintf("Total features analyzed: %d", nrow(feature_matrix)))
logmessage(sprintf("Total samples: %d", ncol(feature_matrix)))
logmessage(sprintf("Optimal clusters identified: %d", optimal_clusters))
logmessage(sprintf("Clustering method: Ward.D2 (Euclidean distance)")
logmessage(sprintf("Average silhouette score: %.4f", avg_silhouette))
logmessage("")
logmessage("Output files:")
logmessage(sprintf(" - Heatmap: %s", file.path(OUTPUT_DIR, "hierarchical_clustering_heatmap.pdf")))
logmessage(sprintf(" - Dendrogram: %s", file.path(OUTPUT_DIR, "dendrogram_colored.pdf")))
logmessage(sprintf(" - Assignments: %s", file.path(OUTPUT_DIR, "cluster_assignments.csv")))
logmessage("="*70)
