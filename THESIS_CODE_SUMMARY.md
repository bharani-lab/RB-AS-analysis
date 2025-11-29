# Thesis Code Repository - Brief Summary

## 📚 Purpose

This GitHub repository contains **all computational code and scripts** used in my PhD thesis:

**"Identification and Analysis of Alternative Transcripts in Retinoblastoma"**

Submitted to Alagappa University, Department of Bioinformatics

Advisor: Dr. D. Bharanidharan

---

## 🎯 Thesis Focus

**Research Question:** What alternative splicing events distinguish retinoblastoma (RB) from normal tissue?

**Approach:** Comprehensive computational analysis using RNA-seq data from 132 RB tumors + 57 fetal controls

**Key Discovery:** Identified 67 RB-specific alternative splicing events via machine learning (Boruta algorithm)

**Validation:** Experimental confirmation through RT-qPCR, minigene assays, and mass spectrometry proteomics

---

## 📂 Code Organization

### Three Thesis Chapters → Computational Pipelines

**Chapter 1: Methods** (RNA-seq data processing)
- Data download and quality control
- Alignment (STAR, HISAT2 benchmarking)
- Transcript quantification (Salmon)
- Alternative splicing detection (SUPPA2, rMATS)

**Chapter 2: Results** (Analysis & discovery)
- Differential alternative splicing (DAS) analysis
- Machine learning feature selection (Boruta - 67 events identified)
- Enrichment analysis (GO, KEGG, Hallmarks)
- Co-expression network analysis (WGCNA)
- Hub gene identification

**Chapter 3: Validation** (Experimental confirmation)
- RT-qPCR validation
- Minigene reporter assays
- Protein structure analysis (AlphaFold2, PyMOL)
- Proteomics validation (LC-MS/MS)
- Splicing factor-target interactions

---

## 🔧 Technologies Used

**RNA-seq Tools:**
- STAR, HISAT2 (alignment)
- Salmon (transcript quantification)
- SUPPA2, rMATS (AS detection)
- fastp (quality control)

**Analysis & Statistics:**
- R (edgeR, Boruta, WGCNA, clusterProfiler)
- Python (pandas, scikit-learn)
- Jupyter notebooks for exploratory analysis

**Visualization & Validation:**
- PyMOL (protein structure)
- Cytoscape (networks)
- ggplot2, ComplexHeatmap (R)
- AlphaFold2 (protein prediction)

---

## 📊 Dataset Summary

| Metric | Value |
|--------|-------|
| **RB Tumor Samples** | 50 |
| **Control Samples (Fetal Retina)** | 17 |
| **Total RNA-seq Samples** | 67 |
| **Sequencing Platform** | Illumina |
| **Reference Genome** | GRCh38 (Ensembl v104) |
| **DAS Events Identified** | 3,408 |
| **RB-Specific Events (Boruta)** | 67 confirmed |
| **Pan-cancer Comparison** | 33 TCGA cancer types |

---

## 📁 Repository Structure

```
RB-AS-analysis/
├── scripts/                 # Analysis scripts (to be added)
│   ├── 01_data_download.sh
│   ├── 02_salmon_quantification.sh
│   ├── 03_boruta_feature_selection.R
│   ├── 04_clustering_heatmap.R
│   └── ...
├── config/                  # Configuration files (to be added)
│   ├── environment.yml
│   ├── parameters.yaml
│   └── sample_metadata.csv
├── docs/                    # Documentation (to be added)
│   ├── methods.md
│   ├── results.md
│   └── validation.md
├── analysis/                # Jupyter notebooks (to be added)
├── results/                 # Output files (to be added)
├── README.md               # Main documentation
└── THESIS_CODE_SUMMARY.md   # This file
```

---

## ⚡ Quick Start

```bash
# Clone repository
git clone https://github.com/bharani-lab/RB-AS-analysis.git
cd RB-AS-analysis

# Setup environment
conda env create -f config/environment.yml
conda activate rb-splicing

# Run analysis pipeline
bash scripts/01_data_download.sh
bash scripts/02_salmon_quantification.sh
Rscript scripts/03_boruta_feature_selection.R
```

---

## 🎓 Citation

If you use this code or analysis:

```bibtex
@phdthesis{aslam2025,
  author = {A. Mohamed Hameed Aslam},
  title = {Identification and Analysis of Alternative Transcripts in Retinoblastoma},
  school = {Alagappa University},
  year = {2025},
  note = {GitHub: github.com/bharani-lab/RB-AS-analysis}
}
```

---

## 📞 Contact

**Student:** A. Mohamed Hameed Aslam (Reg. No: R20162158)

**Advisor:** Dr. D. Bharanidharan  
Department of Bioinformatics  
Aravind Medical Research Foundation, Madurai  

**Email:** bharanidharan@aravind.org

---

## 📄 License

MIT License - Code is freely available for educational and research purposes

---

**Status:** ✅ Repository ready for code upload

**Last Updated:** November 30, 2025
