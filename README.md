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
