# Table of Contents: RB-AS-analysis Scripts

Comprehensive guide to all analysis scripts for retinoblastoma alternative splicing research.

---

## Overview

This directory contains the complete computational pipeline for identifying and validating alternative splicing events in retinoblastoma (RB). Scripts are organized by analysis phase, progressing from raw data processing through to protein-level validation.

**Total Scripts:** 11

**Languages:** Bash (shell), R, Python

**Primary Tools:** HISAT2, Salmon, STAR, rMATS, SUPPA2, Boruta ML, GO/KEGG enrichment, Mass Spectrometry analysis

---

## Script Organization by Analysis Phase

### Phase 1: Data Preprocessing & Quality Control (Scripts 0-3)

Initial RNA-seq data processing and feature selection

#### Script 0: RNA-seq Alignment and Counting
- **File:** `00_rnaseq_alignment_counting.sh`
- **Language:** Bash
- **Description:** Comprehensive RNA-seq alignment pipeline using HISAT2 and feature counting. Aligns raw sequencing reads to reference genome and generates read count matrices.
- **Input:** Raw FASTQ files, reference genome
- **Output:** BAM files, read count matrices, alignment QC metrics
- **Key Tools:** HISAT2, samtools, featureCounts
- **Dependencies:** Reference genome index, GTF annotation

#### Script 1: SRA Data Download
- **File:** `01_sra_data_download.sh`
- **Language:** Bash
- **Description:** Automated download and preprocessing of public RNA-seq data from SRA (Sequence Read Archive). Handles batch downloads and format conversion.
- **Input:** SRA accession numbers
- **Output:** FASTQ files
- **Key Tools:** SRA Toolkit (fastq-dump, prefetch)
- **Use Case:** Acquiring public RB datasets for comparative analysis

#### Script 2: Salmon Transcript Quantification
- **File:** `02_salmon_transcript_quantification.sh`
- **Language:** Bash
- **Description:** Fast and accurate transcript-level quantification using Salmon pseudo-alignment method. Generates transcript abundance estimates (TPM, counts).
- **Input:** FASTQ files, transcript index
- **Output:** Salmon quantification files, TPM/count matrices
- **Key Tools:** Salmon, tximport (R)
- **Advantages:** Fast, lightweight, accurate for isoform abundance
- **SUPPA2 Integration:** Generates PSI calculations for splice variants

#### Script 3: Boruta ML Feature Selection
- **File:** `03_boruta_feature_selection.R`
- **Language:** R
- **Description:** Identifies RB-specific alternative splicing events using the Boruta machine learning algorithm. Performs feature importance ranking on transcript abundance data.
- **Input:** Normalized transcript abundance matrix
- **Output:** Ranked list of important splicing events, feature importance scores
- **Key Tools:** Boruta R package, randomForest
- **Method:** All-relevant feature selection identifying RB-specific transcripts
- **Criteria:** Used to filter events for experimental validation

---

### Phase 2: Exploratory Data Analysis & Visualization (Scripts 4-7)

Comprehensive statistical and graphical analysis of identified splicing events

#### Script 4: Hierarchical Clustering & Heatmap Visualization
- **File:** `04_hierarchical_clustering_heatmap.R`
- **Language:** R
- **Description:** Unsupervised hierarchical clustering of samples and alternative splicing events. Generates publication-quality heatmaps showing sample relationships and event patterns.
- **Input:** Event abundance matrix (PSI values or transcript counts)
- **Output:** Clustered heatmap, dendrogram, sample distance matrix
- **Key Tools:** pheatmap, dendextend, ggplot2
- **Visualization:** Shows tissue-specific and disease-specific clustering patterns

#### Script 5: Violin Plot Comparison Analysis
- **File:** `05_violin_plot_comparison.R`
- **Language:** R
- **Description:** Statistical comparison of splicing event inclusion levels (PSI) between RB and normal retinal tissue. Generates violin plots with significance testing.
- **Input:** PSI values, sample metadata, group information
- **Output:** Violin plots, statistical test results (t-tests, Mann-Whitney U)
- **Key Tools:** ggplot2, tidyverse, ggpubr
- **Statistics:** Includes p-values, effect sizes, confidence intervals

#### Script 6: Box Plot Pan-Cancer Comparison
- **File:** `06_box_plot_pan_cancer.R`
- **Language:** R
- **Description:** Compares RB-specific splicing events across 33 TCGA cancer types. Contextualizes RB events within broader cancer genomics landscape.
- **Input:** RB events, TCGA cancer PSI data (from TCGA sources)
- **Output:** Box plots showing event variation across cancers, cancer-specific signatures
- **Key Tools:** ggplot2, tidyverse, ggpubr
- **Scope:** 33 TCGA cancer types, identifies RB-specific vs. pan-cancer events
- **Interpretation:** Distinguishes RB-specific from general cancer-associated splicing

#### Script 7: GO and KEGG Pathway Enrichment
- **File:** `07_go_kegg_enrichment.R`
- **Language:** R
- **Description:** Functional annotation of genes harboring alternative splicing events. Identifies enriched biological processes, molecular functions, and KEGG pathways.
- **Input:** Gene list (genes with RB-specific AS events), background gene set
- **Output:** Enriched GO terms, KEGG pathways, enrichment plots
- **Key Tools:** clusterProfiler, org.Hs.eg.db, DOSE
- **Analysis:** GO:BP, GO:MF, GO:CC, KEGG pathway enrichment
- **Visualization:** Bubble plots, bar charts, semantic similarity networks

---

### Phase 3: Alternative Splicing Event Detection & Processing (Scripts 8-10)

Core analysis identifying and characterizing alternative splicing events

#### Script 8: Merge Alternative Splicing Event Datasets
- **File:** `08_merge_event_datasets.R`
- **Language:** R
- **Description:** Integrates alternative splicing events detected from multiple methods (rMATS, SUPPA2). Creates unified event matrix for downstream analysis.
- **Input:** rMATS results, SUPPA2 PSI values, event annotation files
- **Output:** Merged event matrix, unified event annotations
- **Key Tools:** tidyverse, data.table
- **Purpose:** Combines complementary detection methods to maximize sensitivity
- **Filtering:** Removes redundant events, applies confidence filters

#### Script 9: Alternative Splicing Event Matching Utility
- **File:** `09_event_matching_utility.R`
- **Language:** R
- **Description:** Utility functions for matching and cross-referencing alternative splicing events across different annotation sources and analysis tools.
- **Input:** Event lists from different sources/formats
- **Output:** Matched event pairs, mapping tables, cross-reference tables
- **Key Functions:** Event coordinate matching, isoform matching, format conversion
- **Use Case:** Standardizing events for comparison across tools and datasets

#### Script 10: Alternative Splicing Detection (rMATS & SUPPA2)
- **File:** `10_rmats_suppa2_alternative_splicing.sh`
- **Language:** Bash
- **Description:** Core alternative splicing detection pipeline combining rMATS and SUPPA2. Identifies exon skipping, intron retention, alternative 5'/3' splice sites, mutually exclusive exons, and other splicing patterns.
- **Input:** BAM files from HISAT2/STAR alignment, GTF annotation
- **Output:** Event tables, PSI values, statistical significance tests
- **Methods:**
  - **rMATS:** Detects 5 major splicing patterns with quantitative metrics
  - **SUPPA2:** Generates per-event PSI for fine-grained analysis
- **Statistics:** FDR-corrected p-values, effect sizes
- **Key Parameters:** Min reads per event, min sample percentage

---

### Phase 4: Protein-Level Validation (Script 13)

Experimental validation at protein level

#### Script 13: Mass Spectrometry Proteomics Analysis
- **File:** `13_mass_spectrometry_proteomics.R`
- **Language:** R
- **Chapter:** Chapter 3: Protein-level validation of alternative splicing
- **Description:** Analyzes mass spectrometry proteomics data to validate alternative splicing events at the protein level. Identifies isoform-specific peptides and quantifies protein abundance changes.
- **Input:** MaxQuant proteomics data (peptide intensity matrix), AS event predictions
- **Output:** Isoform-specific peptide summary, protein coverage analysis, abundance comparisons
- **Key Tools:** MaxQuant PSM data, tidyverse, ggplot2
- **Analysis Type:** Isoform-specific peptide detection and quantification
- **Validation:** Confirms computational AS predictions with protein evidence
- **Integration:** Links transcript-level AS events to protein expression changes

---

## Execution Flow / Pipeline Workflow

```
Phase 1: Data Preprocessing
├── Script 1: Download SRA data → FASTQ files
├── Script 0: Align RNA-seq → BAM files + counts
└── Script 2: Quantify transcripts → TPM/abundance
           ↓
Phase 3: Splicing Event Detection
└── Script 10: rMATS/SUPPA2 → Raw splicing events + PSI
├── Script 8: Merge events → Unified matrix
└── Script 9: Match events → Cross-reference
                     ↓
Phase 2: Feature Selection & Filtering
└── Script 3: Boruta ML → Important splicing events
 ↓
Phase 4: Statistical Analysis & Visualization  
├── Script 4: Hierarchical clustering → Heatmaps
├── Script 5: Violin plots → RB vs normal comparison
├── Script 6: Box plots → Pan-cancer context
└── Script 7: GO/KEGG enrichment → Functional annotation
           ↓
Phase 5: Protein-Level Validation
└── Script 13: Mass spectrometry → Protein isoform validation
```

---

## Input/Output Data Requirements

### Required Input Data
- **Raw RNA-seq:** FASTQ format files
- **Reference genome:** FASTA + GTF annotation (e.g., hg38/GRCh38)
- **Sample metadata:** CSV with sample IDs, tissue type, condition
- **Mass spectrometry data:** MaxQuant output (peptides.txt, Protein groups.txt)

### Output Locations
- **Alignment results:** `results/bam_files/`
- **Quantification:** `results/counts/`, `results/abundance/`
- **Splicing events:** `results/splicing_events/`
- **Plots/Visualizations:** `results/figures/`
- **Statistical results:** `results/statistics/`

---

## Software Dependencies

### Bioinformatics Tools
- **Alignment:** HISAT2, STAR
- **Quantification:** Salmon, featureCounts
- **Splicing Detection:** rMATS, SUPPA2
- **Data Retrieval:** SRA Toolkit
- **File Processing:** samtools, bedtools

### R Packages (Key)
- Data manipulation: tidyverse, data.table
- Visualization: ggplot2, ggpubr, pheatmap
- Machine learning: randomForest, Boruta
- Functional analysis: clusterProfiler, org.Hs.eg.db, DOSE
- Statistics: limma, edgeR, DESeq2 (optional)

### Python Tools
- PyMOL (for protein structure visualization, if used)

---

## Running the Pipeline

### Sequential Execution
```bash
# Phase 1: Data Preprocessing
bash scripts/01_sra_data_download.sh
bash scripts/00_rnaseq_alignment_counting.sh
bash scripts/02_salmon_transcript_quantification.sh

# Phase 2: Feature Selection
Rscript scripts/03_boruta_feature_selection.R

# Phase 3: Splicing Detection
bash scripts/10_rmats_suppa2_alternative_splicing.sh
Rscript scripts/08_merge_event_datasets.R
Rscript scripts/09_event_matching_utility.R

# Phase 4: Analysis & Visualization
Rscript scripts/04_hierarchical_clustering_heatmap.R
Rscript scripts/05_violin_plot_comparison.R
Rscript scripts/06_box_plot_pan_cancer.R
Rscript scripts/07_go_kegg_enrichment.R

# Phase 5: Validation
Rscript scripts/13_mass_spectrometry_proteomics.R
```

---

## Script Statistics

| Phase | Script Count | Languages | Focus |
|-------|---------|-----------|-------|
| Data Preprocessing (0-3) | 4 | Bash, R | Raw data to abundance |
| Analysis & Visualization (4-7) | 4 | R | Statistical analysis, plots |
| Splicing Detection (8-10) | 3 | R, Bash | Core event identification |
| Protein Validation (13) | 1 | R | Proteomics integration |
| **Total** | **11** | **Bash, R, Python** | **Complete pipeline** |

---

## Contact & Citation

**Supervisor:** Dr. D. Bharanidharan, AMRF Lab, Department of Bioinformatics

**Thesis:** Identification and Analysis of Alternative Transcripts in Retinoblastoma Progression

**Author:** A. Mohamed Hameed Aslam (Reg. No: R20162158)

**University:** Allagappa University

---


- 
---

## List of Publications from the Thesis

### Published

- **Assane Rachidou, M. H. A.**, Vanniarajan, A., Kim, U., & Devarajan, B. (2025). Analysing Differential Alternative Splicing Events and Their Impact on Retinoblastoma Progression Using RNA-seq Metadata. *Asian Pacific Journal of Cancer Prevention*, 26(5), 1781–1792. https://doi.org/10.31557/APJCP.2025.26.5.1781

### Under Review

- **Assane Rachidou, M. H. A.**, Sethu, N., Vanniarajan, A., Kim, U., & Devarajan, B. Identification of Differential Alternative Splicing Events Specific to Retinoblastoma Associated with Tumour Progression. *(Under review)*
