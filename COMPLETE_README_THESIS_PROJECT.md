# Computational Analysis of Alternative Splicing in Retinoblastoma

**PhD Thesis Project - Comprehensive RNA-seq & Bioinformatics Pipeline**

🧬 **Research Focus:** Identification and validation of alternative splicing events specific to retinoblastoma using computational and experimental approaches

---

## 📋 Project Overview

### Research Objective

This thesis investigates **alternative splicing (AS) events in retinoblastoma (RB)** progression by:

1. **Characterizing AS landscape** in 50 RB tumors + 17 fetal controls using RNA-seq
2. **Identifying RB-specific DAS events** via machine learning (Boruta feature selection)
3. **Validating findings** through experimental and proteomic analyses
4. **Functional characterization** of splicing factors and target genes

### Key Findings

✅ **67 confirmed AS events** distinguishing RB from control samples
✅ **SF target identification** (ILF2, HNRNPA1, etc.)
✅ **Pathway enrichment** linking AS to cancer hallmarks
✅ **Experimental validation** via RT-qPCR, minigene assays, mass spectrometry

---

## 🗂️ Repository Structure

```
RB-AS-analysis/
├── GITHUB_UPLOAD_GUIDE.md              # Step-by-step upload instructions
├── COMPLETE_README_THESIS_PROJECT.md   # This file
├── LICENSE                             # MIT License
├── .gitignore                          # Python exclusions
│
├── docs/                               # Documentation (for GitHub Pages)
│   ├── index.md                        # Main documentation
│   ├── methods.md                      # Methods details
│   ├── results.md                      # Results summary
│   └── validation.md                   # Validation protocols
│
├── scripts/                            # Analysis scripts (to be added)
│   ├── 01_data_download.sh             # SRA download
│   ├── 02_salmon_quantification.sh     # Transcript quantification
│   ├── 03_boruta_feature_selection.R   # ML feature selection
│   ├── 04_clustering_heatmap.R         # Visualization
│   ├── 05_enrichment_analysis.R        # Pathway analysis
│   ├── 06_network_analysis.R           # Network construction
│   ├── 07_validation_analysis.R        # Experimental validation
│   ├── 08_protein_analysis.R           # Proteomics integration
│   └── 09_visualization.R              # Publication figures
│
├── config/                             # Configuration files (to be added)
│   ├── environment.yml                 # Conda dependencies
│   ├── parameters.yaml                 # Analysis parameters
│   └── sample_metadata.csv             # Sample information
│
├── analysis/                           # Jupyter notebooks (to be added)
│   ├── 01_exploratory_analysis.ipynb
│   ├── 02_das_discovery.ipynb
│   ├── 03_feature_importance.ipynb
│   └── 04_functional_networks.ipynb
│
└── results/                            # Output files (to be added)
    ├── figures/                        # Publication figures
    ├── tables/                         # Summary statistics
    └── README.md                       # Results guide
```

---

## 🔬 Thesis Structure (3 Chapters)

### Chapter 1: Methods & Datasets
**"Computational Identification of Alternative Splicing Events in Retinoblastoma"**

- RNA-seq data collection (67 samples: 50 RB, 17 controls)
- Quality control & preprocessing
- Alignment & quantification (STAR, Salmon)
- AS event detection (SUPPA2, rMATS)
- References: `scripts/01-02`, `docs/methods.md`

### Chapter 2: Results & Analysis
**"RB-Specific Alternative Splicing Signature Identification"**

- DAS event landscape characterization
- Machine learning feature selection (Boruta)
- 67 confirmed RB-specific events
- Enrichment analysis (GO, KEGG, Hallmarks)
- Hub gene identification (WGCNA)
- Network analysis & visualization
- References: `scripts/03-06`, `docs/results.md`

### Chapter 3: Validation & Functional Analysis
**"Experimental Validation & Functional Impact of AS Events"**

- RT-qPCR validation of selected events
- Minigene reporter assays
- Proteomics validation (LC-MS/MS)
- SF-target interactions
- Clinical implications
- References: `scripts/07-09`, `docs/validation.md`

---

## 📊 Dataset Summary

| Component | Value | Notes |
|-----------|-------|-------|
| RB tumors | 50 | SRA accessions available |
| Control samples | 17 | Fetal retinal tissue |
| Total samples | 67 | RNA-seq (Illumina) |
| Reference genome | GRCh38 v104 | Ensembl annotation |
| DAS events identified | 67 | After filtering & ML |
| Confirmed events | 67 | Boruta importance score |

---

## 🛠️ Required Tools & Dependencies

### Bioinformatics Tools

```yaml
# RNA-seq processing
- STAR: 2.7.9a (alignment)
- HISAT2: 2.2.1 (alternative aligner)
- Salmon: 1.6.0 (transcript quantification)
- fastp: Quality control

# Alternative splicing
- SUPPA2: DAS detection
- rMATS: AS event annotation

# Analysis & visualization
- R 4.1+ (base statistical environment)
- Python 3.9+ (bioinformatics scripting)
- Jupyter Lab (interactive analysis)

# R packages
- edgeR, DESeq2 (differential expression)
- Boruta, randomForest (ML feature selection)
- WGCNA (co-expression networks)
- clusterProfiler (enrichment analysis)
- igraph, tidygraph (network analysis)
- ggplot2, ComplexHeatmap (visualization)

# Python packages
- pandas, numpy, scipy (data analysis)
- scikit-learn (machine learning)
- networkx, igraph (network analysis)
- matplotlib, seaborn, plotly (visualization)
```

### Quick Setup with Conda

```bash
# Create environment
conda env create -f config/environment.yml
conda activate rb-splicing

# Verify installation
R --version
python --version
STAR --version
salmon --version
```

---

## 🚀 Quick Start Guide

### 1. **Clone Repository**

```bash
git clone https://github.com/bharani-lab/RB-AS-analysis.git
cd RB-AS-analysis
```

### 2. **Setup Environment**

```bash
conda env create -f config/environment.yml
conda activate rb-splicing
```

### 3. **Review Documentation**

- Start with `GITHUB_UPLOAD_GUIDE.md` for upload workflow
- See `docs/methods.md` for detailed methodology
- Check `docs/results.md` for findings
- Review `docs/validation.md` for experimental protocols

### 4. **Run Analysis Scripts** (after scripts are uploaded)

```bash
# Data preprocessing
bash scripts/01_data_download.sh

# Quantification
bash scripts/02_salmon_quantification.sh

# Feature selection
Rscript scripts/03_boruta_feature_selection.R

# Visualization
Rscript scripts/04_clustering_heatmap.R
```

### 5. **Explore Notebooks**

```bash
jupyter lab analysis/01_exploratory_analysis.ipynb
```

---

## 📁 Files Available for Download

From your Perplexity thread, these comprehensive guides are ready:

✅ **scripts_with_comments_part1.md** (50+ KB)
   - Script 1: `data_download.sh` with 200+ comment lines
   - Script 2: `salmon_transcript_quant.sh` with 170+ comment lines

✅ **scripts_with_comments_part2.md** (50+ KB)
   - Script 3: `boruta_feature_selection.R` with 280+ comment lines
   - Script 4: `hierarchical_clustering_heatmap.R` with 240+ comment lines

✅ **github_upload_guide_part3.md** (40+ KB)
   - Scripts 5-9 complete templates
   - Master upload workflow (PHASE 1-6)

✅ **FINAL_COMPLETE_DELIVERY_SUMMARY.md** (20+ KB)
   - Complete overview & reference

✅ **github_pages_complete_guide.md** (50+ KB)
   - GitHub Pages setup for thesis publication

---

## 📌 File Upload Status

| File | Status | Location | Notes |
|------|--------|----------|-------|
| GITHUB_UPLOAD_GUIDE.md | ✅ Uploaded | Root | Phase-by-phase instructions |
| COMPLETE_README_THESIS_PROJECT.md | ✅ Uploaded | Root | This comprehensive overview |
| config/environment.yml | ⏳ Ready | config/ | Download from Perplexity |
| config/parameters.yaml | ⏳ Ready | config/ | Download from Perplexity |
| scripts/*.sh, *.R | ⏳ Ready | scripts/ | Extract from Perplexity files |
| docs/*.md | ⏳ Ready | docs/ | Use templates from guides |

---

## 🔗 GitHub Links

- **Repository:** https://github.com/bharani-lab/RB-AS-analysis
- **GitHub Pages:** https://bharani-lab.github.io/RB-AS-analysis (when enabled)
- **Perplexity Guides:** https://www.perplexity.ai/search/i-want-to-upload-all-my-codes-GkXoO3ILS12O8Zz2AnonUw

---

## 📖 Citation

If you use this code or analysis, please cite:

```bibtex
@thesis{thesis2025,
  author = {[Your Name]},
  title = {Computational Analysis of Alternative Splicing in Retinoblastoma},
  school = {Allagappa University},
  year = {2025},
  note = {GitHub: github.com/bharani-lab/RB-AS-analysis}
}
```

---

## 👥 Contact & Support

**Advisor:** Dr. Bharanidharan  
**Department:** Department of Bioinformatics  
**Institution:** Allagappa University / AMRF Lab

**Questions or Issues?** 
- Check `GITHUB_UPLOAD_GUIDE.md` for troubleshooting
- Review `docs/faq.md` (when added) for common questions
- See specific script documentation for usage help

---

## 📜 License

This project is licensed under the **MIT License** - see LICENSE file for details.

This allows open science sharing while protecting your intellectual property.

---

## ✅ Project Status

| Component | Status | Notes |
|-----------|--------|-------|
| GitHub Repository | ✅ Active | Public repository created |
| Upload Guide | ✅ Complete | GITHUB_UPLOAD_GUIDE.md added |
| Core Documentation | ✅ Complete | README & project overview |
| Scripts (ready to upload) | ✅ Prepared | 9 fully commented scripts |
| Config Files | ✅ Prepared | environment.yml, parameters.yaml |
| GitHub Pages | ⏳ Ready | Can enable in Settings → Pages |
| Thesis Publication | ⏳ Pending | Follow upload guide workflow |

---

## 🎓 Next Steps

1. **Follow GITHUB_UPLOAD_GUIDE.md** for detailed upload instructions
2. **Extract scripts** from Perplexity files (scripts_with_comments_part*.md)
3. **Create directory structure** (scripts/, config/, docs/, analysis/)
4. **Upload scripts** following Phase 1-6 workflow
5. **Configure GitHub Pages** for public thesis documentation
6. **Archive on Zenodo** for permanent DOI and citation

---

**Last Updated:** November 30, 2025

**Repository Status:** ✅ Ready for thesis code upload

**Next Action:** Download guides from Perplexity and follow GITHUB_UPLOAD_GUIDE.md
