#!/usr/bin/env Rscript

# Script 7: GO & KEGG Pathway Enrichment Analysis
# Purpose: Perform Gene Ontology and KEGG pathway enrichment on RB alternative splicing events
# Author: Bharani Lab
# Date: 2024

# Load required libraries
require(clusterProfiler)
require(org.Hs.eg.db)
require(enrichplot)
require(data.table)
require(dplyr)
require(ggplot2)

# Configuration
INPUT_FILE <- "data/confirmed_rb_features.RDS"
GENE_MAPPING <- "data/gene_id_mapping.csv"
OUTPUT_DIR <- "results/07_enrichment_analysis"

if (!dir.exists(OUTPUT_DIR)) {
  dir.create(OUTPUT_DIR, recursive = TRUE)
}

# ============================================================================
# 1. LOAD DATA AND PREPARE GENE IDS
# ============================================================================
logmessage("Loading RB feature data...")
confirmed_features <- readRDS(INPUT_FILE)

logmessage("Loading gene ID mapping...")
gene_mapping <- fread(GENE_MAPPING)

# Extract feature names and map to Entrez IDs
feature_names <- confirmed_features$feature_name
feature_to_gene <- merge(
  data.table(feature = feature_names),
  gene_mapping,
  by.x = "feature", by.y = "feature_name",
  all.x = TRUE
)

entrez_ids <- feature_to_gene$entrez_id
entrez_ids <- na.omit(entrez_ids)

logmessage(sprintf("Mapped %d features to %d unique genes", length(feature_names), length(unique(entrez_ids))))

# ============================================================================
# 2. GO ENRICHMENT ANALYSIS
# ============================================================================
logmessage("Performing GO enrichment analysis...")

# Biological Process
go_bp <- enrichGO(
  gene = entrez_ids,
  OrgDb = org.Hs.eg.db,
  ont = "BP",
  pAdjustMethod = "BH",
  pvalueCutoff = 0.05,
  readable = TRUE
)

logmessage(sprintf("GO-BP: %d significantly enriched terms", nrow(go_bp@result)))

# Molecular Function
go_mf <- enrichGO(
  gene = entrez_ids,
  OrgDb = org.Hs.eg.db,
  ont = "MF",
  pAdjustMethod = "BH",
  pvalueCutoff = 0.05,
  readable = TRUE
)

logmessage(sprintf("GO-MF: %d significantly enriched terms", nrow(go_mf@result)))

# Cellular Component
go_cc <- enrichGO(
  gene = entrez_ids,
  OrgDb = org.Hs.eg.db,
  ont = "CC",
  pAdjustMethod = "BH",
  pvalueCutoff = 0.05,
  readable = TRUE
)

logmessage(sprintf("GO-CC: %d significantly enriched terms", nrow(go_cc@result)))

# ============================================================================
# 3. KEGG PATHWAY ANALYSIS
# ============================================================================
logmessage("Performing KEGG pathway enrichment analysis...")

kegg_result <- enrichKEGG(
  gene = entrez_ids,
  organism = "hsa",
  pAdjustMethod = "BH",
  pvalueCutoff = 0.05
)

logmessage(sprintf("KEGG: %d significantly enriched pathways", nrow(kegg_result@result)))

if (nrow(kegg_result@result) > 0) {
  kegg_result <- setReadable(kegg_result, OrgDb = org.Hs.eg.db, keyType = "ENTREZID")
}

# ============================================================================
# 4. EXPORT ENRICHMENT RESULTS
# ============================================================================
logmessage("Exporting enrichment results...")

# Export GO results
fwrite(as.data.table(go_bp@result), file.path(OUTPUT_DIR, "GO_biological_process.csv"))
fwrite(as.data.table(go_mf@result), file.path(OUTPUT_DIR, "GO_molecular_function.csv"))
fwrite(as.data.table(go_cc@result), file.path(OUTPUT_DIR, "GO_cellular_component.csv"))

# Export KEGG results
if (nrow(kegg_result@result) > 0) {
  fwrite(as.data.table(kegg_result@result), file.path(OUTPUT_DIR, "KEGG_pathways.csv"))
}

logmessage("Enrichment results exported")

# ============================================================================
# 5. GENERATE VISUALIZATIONS
# ============================================================================
logmessage("Generating enrichment visualizations...")

# GO Bar plots
if (nrow(go_bp@result) > 0) {
  p_bp <- barplot(go_bp, showCategory = 15, title = "GO Biological Process Enrichment")
  ggsave(file.path(OUTPUT_DIR, "GO_BP_barplot.pdf"), p_bp, width = 12, height = 8)
}

if (nrow(go_mf@result) > 0) {
  p_mf <- barplot(go_mf, showCategory = 15, title = "GO Molecular Function Enrichment")
  ggsave(file.path(OUTPUT_DIR, "GO_MF_barplot.pdf"), p_mf, width = 12, height = 8)
}

# KEGG plot
if (nrow(kegg_result@result) > 0) {
  p_kegg <- barplot(kegg_result, showCategory = 15, title = "KEGG Pathway Enrichment")
  ggsave(file.path(OUTPUT_DIR, "KEGG_barplot.pdf"), p_kegg, width = 12, height = 8)
  
  # Dotplot
  p_kegg_dot <- dotplot(kegg_result, showCategory = 15)
  ggsave(file.path(OUTPUT_DIR, "KEGG_dotplot.pdf"), p_kegg_dot, width = 10, height = 8)
}

# Network plot (if sufficient terms)
if (nrow(go_bp@result) >= 10) {
  p_net <- emapplot(go_bp, showCategory = 30)
  ggsave(file.path(OUTPUT_DIR, "GO_BP_network.pdf"), p_net, width = 14, height = 12)
}

logmessage("Visualizations saved")

# ============================================================================
# 7. SUMMARY STATISTICS
# ============================================================================
logmessage("Creating enrichment summary...")

enrichment_summary <- data.table(
  Analysis = c("GO-BP", "GO-MF", "GO-CC", "KEGG"),
  Enriched_Terms = c(nrow(go_bp@result), nrow(go_mf@result), nrow(go_cc@result), nrow(kegg_result@result)),
  Top_P_Value = c(
    ifelse(nrow(go_bp@result) > 0, min(go_bp@result$p.adjust), NA),
    ifelse(nrow(go_mf@result) > 0, min(go_mf@result$p.adjust), NA),
    ifelse(nrow(go_cc@result) > 0, min(go_cc@result$p.adjust), NA),
    ifelse(nrow(kegg_result@result) > 0, min(kegg_result@result$p.adjust), NA)
  )
)

fwrite(enrichment_summary, file.path(OUTPUT_DIR, "enrichment_summary.csv"))

# ============================================================================
# OUTPUT SUMMARY
# ============================================================================
logmessage("")
logmessage("="*70)
logmessage("GO & KEGG ENRICHMENT ANALYSIS COMPLETE")
logmessage("="*70)
logmessage(sprintf("Genes analyzed: %d", length(unique(entrez_ids))))
logmessage(sprintf("GO-BP enriched terms: %d", nrow(go_bp@result)))
logmessage(sprintf("GO-MF enriched terms: %d", nrow(go_mf@result)))
logmessage(sprintf("GO-CC enriched terms: %d", nrow(go_cc@result)))
logmessage(sprintf("KEGG enriched pathways: %d", nrow(kegg_result@result)))
logmessage("")
logmessage("Output files generated in:", OUTPUT_DIR)
logmessage("="*70)
