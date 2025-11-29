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

### Problem Statement
- Retinoblastoma (~1 in 20,000 births) - most common childhood eye malignancy
- 4-15% of RB1 mutations occur at splice sites, suggesting abnormal splicing drives progression
- Previous studies limited to single datasets; comprehensive pan-cancer comparison was missing
- Critical gap: RB splicing event specificity versus other cancers unknown

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
 - **Critical p53 inactivation pathway:** MDM2/MDM4 splicing dysregulation (RPS15, RPL37, CDK5RAP3)
  
### Chapter 3: Experimental Validation
- **RT-qPCR validation:** 5 RB tumors + 4 retinal controls (RNA integrity 7.2-8.9)
- **Statistically significant events:**
  - CCNB1 exon-skip: p=0.01, fold change 6.02
  - ENO2 exon-skip: p=0.01, fold change 4.91
  - CDK5RAP3 intron retention: p=0.03, fold change 4.5
- **Proteomics evidence:** CDK5RAP3 20.72% sequence coverage (2 unique peptides, 7 total PSMs)

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
- 
### Additional Structural Impact Analysis

Other splicing-affected proteins exhibited distinct domain-level perturbations consistent with functional disruption:

#### RB1 (Pocket/Linker Region Perturbation)

- **Mechanism:** Splice-site-altering variants affecting exons enriched around the A/B pocket and adjacent linker segments distort the cyclin-like pocket architecture critical for E2F binding and chromatin tethering
- **Structural consequence:** Weakened pRB-mediated transcriptional repression and checkpoint control, complementing classical truncating RB1 mutations
- **Impact on RB progression:** Enhanced loss of growth suppression function in RB1-mutant tumours

#### PCNA (DNA Clamp Surface Rewiring)

- **Mechanism:** Alternative splicing events mapping to surface loops of PCNA modify the electrostatic landscape of the DNA clamp and partner-binding interfaces
- **Structural consequence:** Altered affinity for replication and repair factors including pRB-regulated complexes; shifts balance between high-fidelity repair and error-prone DNA synthesis
- **Impact on RB progression:** Promotes genomic instability and enhanced mutation burden in RB cells

#### TFDP1 (E2F Dimerization/DNA-Recognition Domain)

- **Mechanism:** TFDP1 splice variants impact the dimerization or DNA-recognition region of the DP1 subunit, perturbing E2F–DP complex assembly on target promoters
- **Structural consequence:** Reconfigures E2F target gene programmes controlling cell cycle, apoptosis, and DNA repair
- **Impact on RB progression:** Reinforces proliferative output of RB1 pathway disruption and altered E2F-dependent transcriptome

#### Splicing Factor Nodes (HNRNPA1, ILF2, SR Proteins)

- **Mechanism:** Differential splicing of core RNA-binding proteins and SR/hnRNP factors alters low-complexity and RS-rich regions mediating RNA recognition and multivalent protein–protein interactions
- **Structural consequence:** Changes in these intrinsically disordered or modular domains propagate global splicing pattern shifts, indirectly reshaping domain composition and interactomes of numerous tumour-relevant proteins
- **Impact on RB progression:** Establishes sustained splicing dysregulation network affecting multiple cell cycle checkpoint nodes simultaneously

#### Integrative Model: Splicing-Driven Secondary Hit Hypothesis

These structural alterations establish that alternative splicing creates protein isoforms with defined domain disruptions affecting:
- Cell cycle regulation (CCNB1 NES disruption, CDK5RAP3 checkpoint loss)
- Metabolic enzyme activity (ENO2 active site alteration)
- Splicing factor function (LUC7L RS-domain loss)
- Transcriptional control (TFDP1 E2F interaction, RB1 pocket distortion)
- DNA replication fidelity (PCNA clamp restructuring)

The convergence of splicing-driven alterations across multiple cell cycle control nodes creates a selective advantage for tumour cells through coordinated dysregulation, representing a complementary mechanism to RB1 mutations and defining alternative splicing as both a critical vulnerability and therapeutic opportunity.


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

