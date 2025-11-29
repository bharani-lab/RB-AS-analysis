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
- 
- **Critical p53 inactivation pathway:** MDM2/MDM4 splicing dysregulation (RPS15, RPL37, CDK5RAP3)
- 
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

## Key Findings Summary

### Problem Statement
- Retinoblastoma (~1 in 20,000 births) - most common childhood eye malignancy
- 4-15% of RB1 mutations occur at splice sites, suggesting abnormal splicing drives progression
- Previous studies limited to single datasets; comprehensive pan-cancer comparison was missing
- Critical gap: RB splicing event specificity versus other cancers unknown

### Protein Structure Impacts & Functional Consequences

#### CCNB1 (Cell Cycle Regulator - NES Domain Disruption)
- **Validation:** RT-qPCR p=0.01, fold change 6.02
- **Structural alterations:** RMSD 1.87 A, TM-score 0.8, surface area change 230 sq A
- **Domain disrupted:** Nuclear Export Signal (NES) - critical for subcellular localization
- **Mechanism:** Exon skipping replaces α-helix with loop, junction displacement 4.2 A
- **Consequence:** Aberrant cyclin B1 accumulation → uncontrolled mitotic progression

#### ENO2 (Metabolic Enzyme - Active Site Alteration)
- **Validation:** RT-qPCR p=0.01, fold change 4.91
- **Structural alterations:** RMSD 1.23 A, TM-score 0.88, surface area change 110 sq A
- **Domain disrupted:** Active site and substrate binding pocket
- **Mechanism:** Loss of β-sheet repositions substrate binding, junction displacement 2.1 A
- **Consequence:** Loss of glycolytic function → metabolic reprogramming toward proliferation

#### LUC7L (RNA Processing Factor - RS-Domain Loss)
- **Validation:** RT-qPCR p=0.9 (negative - tissue-specific, not disease-specific)
- **Structural alterations:** RMSD 2.22 A, TM-score 0.48, surface area change 305 sq A
- **Domain disrupted:** Complete RS-domain deletion - largest surface change observed
- **Mechanism:** Junction displacement 6.3 A, removes protein-RNA interaction capability
- **Consequence:** Loss of spliceosome assembly → dominant-negative regulatory effects

#### CDK5RAP3 (Checkpoint Control - LXXLL Motif Disruption)
- **Validation:** RT-qPCR p=0.03, fold change 4.5 | Proteomics: 20.72% sequence coverage
- **Structural alterations:** RMSD 5.18 A, TM-score 0.42 (most dramatic), surface area 160 sq A
- **Domain disrupted:** LXXLL motif - nuclear receptor/regulatory protein interaction
- **Mechanism:** Intron retention disrupts LXXLL helix, junction displacement 3.8 A
- **Consequence:** Elevated CDK5 activity → aberrant progression, reduced checkpoint control

### Multi-Layered Tumor Suppressor Inactivation
- **Layer 1:** Classical RB1 loss disrupts E2F pathway cell cycle control
- **Layer 2:** MDM2/MDM4 splicing dysregulation impairs p53 surveillance
  - RPS15 (alt 5' splice site), RPL37 (skipped exon), CDK5RAP3 (intron retention)
- **Outcome:** Redundant inactivation resistant to single-pathway interventions

### Disease-Specific Molecular Signatures
- **LUC7L:** PSI -0.18 in RB vs +0.12 glioblastoma, +0.09 HCC
- **CDK5RAP3:** PSI +0.22 in RB vs -0.15 lung cancer, -0.11 breast cancer
- **Clinical implication:** Requires precision medicine tailored to RB's unique profile

### Convergence on Cell Cycle Control
- Alternative splicing creates isoforms with specific domain disruptions:
  - CCNB1 NES dysfunction (cell cycle progression)
  - ENO2 active site loss (metabolic activity)
  - LUC7L spliceosome interaction loss (RNA processing)
  - CDK5RAP3 interaction loss (checkpoint control)
- **Selective advantage:** Coordinated dysregulation = complementary mechanism to RB1 mutations

### Therapeutic & Diagnostic Opportunities
- **Alternative splicing as vulnerability:** Multiple dysregulated nodes create targetable pathways
- **Diagnostic biomarkers:** RB-specific events enable early detection and patient stratification
- **Structure-based therapeutics:** Domain disruptions provide targets for isoform-specific modulation
- **Precision stratification:** 3 molecular subtypes support personalized treatment planning
