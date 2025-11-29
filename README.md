# RB-AS-analysis
Computational analysis pipeline and code for identifying and validating alternative splicing events in retinoblastoma. RNA-seq analysis, machine learning feature selection, protein structure analysis, and functional validation.

## About This Research

**Thesis:** Identification and Analysis of Alternative Transcripts in Retinoblastoma Progression

**Author:** A. Mohamed Hameed Aslam (Reg. No: R20162158)

**Advisor:** Dr. D. Bharanidharan, AMRF Lab, Department of Bioinformatics, Aravind Medical Research Foundation

**University:** Alagappa University

## Research Overview

This repository contains computational code and analysis pipelines developed for PhD thesis research on alternative splicing in retinoblastoma (RB), the most common eye cancer in children. The research identifies and validates RB-specific alternative splicing events through comprehensive transcriptomic analysis and machine learning approaches.

## Key Findings

- **Dataset:** 67 RNA-seq samples (50 RB tumors + 17 fetal retinal controls)
- **Differential Alternative Splicing Events:** 6,136 DAS events identified, with 67 confirmed RB-specific events using Boruta machine learning algorithm
- **Event Types:** Primarily mutually exclusive exons (41.6%) and exon skipping (35.8%)
- **Validation:** Experimental confirmation through RT-qPCR, mass spectrometry, and protein structure analysis

- ### Chapter 1: Comprehensive Meta-Analysis
- **6,136 significant differential alternative splicing (DAS) events** identified across 1,262 distinct genes
- **Mutually Exclusive Exons (MXE):** 41.6% (2,553 events) - predominant form of dysregulation
- **Skipped Exons (SE):** 35.8% (2,195 events) - broader gene-level impact affecting 734 genes
- **Retained Introns (RI):** 8.4% (515 events)
- **Alternative splice sites:** 8.1% (A3SS) and 6.1% (A5SS) respectively
- **Optimal tool performance:** STAR (98% precision, 92% junction detection), HISAT2 (92% junction accuracy), Hera-EM (r=0.93)
- **Enriched pathways:** RNA splicing regulation, metabolic reprogramming, cell cycle control, oxidative phosphorylation
- **Hub genes:** TFDP1, PCNA, CCNB1 (strongest associations with RB progression)
- **Master regulators:** ILF2 and HNRNPA1 coordinate alternative splicing machinery

### Chapter 2: Machine Learning-Based RB-Specific Identification
- **92 confirmed RB-specific alternative splicing events** identified using Boruta algorithm
- **Refined from 3,408 DAS events** via random forest classification against 9,437 TCGA tumors (33 cancer types)
- **Exon skipping dominance:** 55% of RB-specific events - disease-specific targeting of exon recognition
- **3 molecular subtypes** identified based on splicing profiles - indicates disease heterogeneity
- **Critical p53 inactivation pathway:** MDM2/MDM4 splicing dysregulation (RPS15, RPL37, CDK5RAP3)
- **Disease specificity:** LUC7L PSI -0.18 in RB vs +0.12 in glioblastoma; CDK5RAP3 PSI +0.22 in RB vs -0.15 in lung cancer

### Chapter 3: Experimental Validation
- **RT-qPCR validation:** 5 RB tumors + 4 retinal controls (RNA integrity 7.2-8.9)
- **Statistically significant events:**
  - CCNB1 exon-skip: p=0.01, fold change 6.02
  - ENO2 exon-skip: p=0.01, fold change 4.91
  - CDK5RAP3 intron retention: p=0.03, fold change 4.5
- **Proteomics evidence:** CDK5RAP3 20.72% sequence coverage (2 unique peptides, 7 total PSMs)
- **Structural impact:** CCNB1 NES domain disruption (RMSD 1.87, TM-score 0.8) - impairs nuclear-cytoplasmic shuttling

## Thesis Structure

1. **Chapter 1:** Comprehensive analysis of alternative splicing in RB using meta-analysis of multi-project RNA-seq data
2. **Chapter 2:** Machine learning-based identification of RB-specific alternative splicing events
3. **Chapter 3:** Experimental validation through RT-qPCR, proteomics, and structural characterization

## Technologies Used

- RNA-seq alignment: STAR, HISAT2
- Transcript quantification: Salmon, featureCounts
- Splicing analysis: rMATS, SUPPA2
- Machine learning: Boruta algorithm, Random Forest
- Network analysis: WGCNA
- Proteomics validation: Mass spectrometry (LC-MS/MS)
- Protein structure: AlphaFold2, PyMOL
- Functional enrichment: Gene Ontology, KEGG, MSigDB

## Repository Purpose

This repository serves to archive and share the computational code and scripts used in all three thesis chapters, enabling reproducibility of the analysis and supporting academic publication and committee review.
