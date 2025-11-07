# BCL2A1-Based Tri-Marker Signature for ICB Response in Lung Adenocarcinoma

[![DOI](https://zenodo.org/badge/DOI/10.5281/zenodo.XXXXXXX.svg)](https://doi.org/10.5281/zenodo.XXXXXXX)
[![License: MIT](https://img.shields.io/badge/License-MIT-yellow.svg)](https://opensource.org/licenses/MIT)
[![R Version](https://img.shields.io/badge/R-%E2%89%A54.3.0-blue.svg)](https://www.r-project.org/)
[![Python Version](https://img.shields.io/badge/Python-%E2%89%A53.9-blue.svg)](https://www.python.org/)

## Overview

This repository contains all analysis code and documentation for the manuscript:

**"Integrative Computational Analysis Identifies BCL2A1 as a CD8⁺ T-cell Survival Marker Associated with Immunotherapy Response in Lung Adenocarcinoma"**

**Authors:** Hoang Minh Quan Pham¹'²'³, Po-Hao Feng⁴'⁵, Chia-Ling Chen⁶, Kang-Yun Lee⁴'⁵'⁷'⁸, Chiou-Feng Lin¹'³'⁸'⁹'¹⁰*

**Affiliations:**
1. International Ph.D. Program in Medicine, College of Medicine, Taipei Medical University
2. Department of Oncology, Faculty of Medicine, Can Tho University of Medicine and Pharmacy
3. Department of Microbiology and Immunology, School of Medicine, Taipei Medical University
4. Division of Pulmonary Medicine, Department of Internal Medicine, Shuang Ho Hospital
5. Division of Pulmonary Medicine, School of Medicine, Taipei Medical University
6. School of Respiratory Therapy, College of Medicine, Taipei Medical University
7. Graduate Institute of Clinical Medicine, College of Medicine, Taipei Medical University
8. TMU Research Center of Thoracic Medicine, Taipei Medical University
9. Graduate Institute of Medical Sciences, College of Medicine, Taipei Medical University
10. Core Laboratory of Immune Monitoring, Office of Research & Development, Taipei Medical University

**Corresponding Author:** Prof. Chiou-Feng Lin (cflin2014@tmu.edu.tw)

**Published in:** Computational and Structural Biotechnology Journal (2025)

---

## Table of Contents

- [Overview](#overview)
- [Key Findings](#key-findings)
- [Citation](#citation)
- [Repository Structure](#repository-structure)
- [Installation](#installation)
- [Data Availability](#data-availability)
- [Quick Start](#quick-start)
- [Module Descriptions](#module-descriptions)
- [Reproducibility](#reproducibility)
- [System Requirements](#system-requirements)
- [Troubleshooting](#troubleshooting)
- [Contributing](#contributing)
- [Contact](#contact)
- [License](#license)
- [Acknowledgments](#acknowledgments)

---

## Key Findings

This study establishes **BCL2A1** as a novel biomarker for immune checkpoint blockade (ICB) response in lung adenocarcinoma through:

✨ **Multi-Cohort Discovery & Validation**
- Discovery: n=60 patients (GSE161537)
- Validation: n=63 patients across 3 independent cohorts
- BCL2A1 associated with improved survival in ICB (HR=0.43) but not in non-ICB settings (HR=1.07)

🧬 **Single-Cell Mechanistic Insights**
- 130,584 CD8⁺ T cells analyzed (GSE176021)
- BCL2A1 enriched in tissue-resident memory (Trm) and proliferating subsets
- Preferential expression in pathological responders

📊 **Superior Predictive Model**
- Tri-marker model (BCL2A1 + CD274 + HOT): AUC=0.826 (discovery), 0.757 (validation)
- Outperforms PD-L1 alone (AUC=0.557)
- Demonstrated clinical utility via decision curve analysis

🔬 **Cross-Platform Validation**
- In silico simulation across 3 diagnostic platforms
- NanoString nCounter: ρ=0.993 (optimal concordance)
- HTG EdgeSeq: ρ=0.986 (cost-effective alternative)
- RT-qPCR: ρ=0.982 (batch-sensitive)

---

## Citation

If you use this code or data in your research, please cite:

```bibtex
@article{pham2025bcl2a1,
  title={Integrative Computational Analysis Identifies BCL2A1 as a CD8+ T-cell Survival Marker Associated with Immunotherapy Response in Lung Adenocarcinoma},
  author={Pham, Hoang Minh Quan and Feng, Po-Hao and Chen, Chia-Ling and Lee, Kang-Yun and Lin, Chiou-Feng},
  journal={Computational and Structural Biotechnology Journal},
  year={2025},
  doi={10.XXXX/XXXXX},
  url={https://github.com/yourusername/BCL2A1-LUAD-ICB-Biomarker}
}
```

---

## Repository Structure

```
BCL2A1-LUAD-ICB-Biomarker/
├── README.md                          # This file
├── LICENSE                            # MIT License
├── CITATION.cff                       # Citation metadata
├── requirements.txt                   # Python dependencies
├── install_packages.R                 # R package installer
├── environment.yml                    # Conda environment
├── .gitignore                         # Git ignore rules
│
├── src/                               # Source code
│   ├── module_0_preliminary/          # Preliminary characterization
│   │   ├── 01_differential_expression.R
│   │   ├── 02_survival_analysis.R
│   │   ├── 03_immune_correlation.R
│   │   ├── 04_pathway_enrichment.R
│   │   ├── 05_meta_analysis.R
│   │   └── 06_deconvolution_epic.R
│   │
│   ├── module_A_discovery/            # Discovery modeling
│   │   ├── 01_data_preprocessing.R
│   │   ├── 02_feature_engineering.R
│   │   ├── 03_model_training.R
│   │   ├── 04_nested_cv.R
│   │   ├── 05_model_evaluation.R
│   │   └── 06_comparative_analysis.R
│   │
│   ├── module_B_validation/           # External validation
│   │   ├── 01_load_validation_cohorts.R
│   │   ├── 02_batch_correction.R
│   │   ├── 03_loco_validation.R
│   │   ├── 04_calibration.R
│   │   └── 05_decision_curve_analysis.R
│   │
│   ├── module_C_scRNAseq/            # Single-cell analysis
│   │   ├── 01_data_loading.py
│   │   ├── 02_quality_control.py
│   │   ├── 03_normalization_scaling.py
│   │   ├── 04_clustering.py
│   │   ├── 05_annotation.py
│   │   ├── 06_bcl2a1_analysis.py
│   │   ├── 07_differential_abundance.py
│   │   └── 08_signature_extraction.py
│   │
│   ├── module_D_deconvolution/       # Bulk deconvolution
│   │   ├── 01_reference_building.R
│   │   ├── 02_epic_deconvolution.R
│   │   ├── 03_enrichment_analysis.R
│   │   ├── 04_ablation_experiments.R
│   │   └── 05_robustness_testing.R
│   │
│   ├── module_E_simulation/          # Cross-platform simulation
│   │   ├── 01_platform_simulation.R
│   │   ├── 02_concordance_analysis.R
│   │   ├── 03_stress_testing.R
│   │   ├── 04_negative_controls.R
│   │   └── 05_platform_comparison.R
│   │
│   └── utils/                         # Utility functions
│       ├── data_loader.R
│       ├── geo_download.R
│       ├── plotting_functions.R
│       ├── statistical_tests.R
│       ├── model_utils.R
│       └── helper_functions.R
│
├── data/                              # Data directory
│   ├── README.md                      # Data description
│   ├── download_data.sh               # Automated download script
│   ├── metadata/                      # Sample metadata
│   │   ├── GSE161537_metadata.csv
│   │   ├── validation_cohorts_metadata.csv
│   │   └── cohort_summary.csv
│   └── raw/                           # Raw data (gitignored)
│
├── results/                           # Analysis outputs
│   ├── figures/                       # Main figures (Fig 1-8)
│   ├── supplementary_figures/         # Supplementary figures
│   ├── tables/                        # Main tables
│   ├── supplementary_tables/          # Supplementary tables (S1-S21)
│   └── models/                        # Saved model objects
│
├── notebooks/                         # Interactive notebooks
│   ├── 00_data_exploration.Rmd
│   ├── 01_module_0_walkthrough.Rmd
│   ├── 02_full_pipeline.Rmd
│   └── 03_reproduce_figures.Rmd
│
└── docs/                              # Documentation
    ├── installation.md
    ├── data_description.md
    ├── analysis_workflow.md
    ├── platform_specifications.md
    └── troubleshooting.md
```

---

## Installation

### Prerequisites

- **R** version ≥ 4.3.0
- **Python** version ≥ 3.9
- **RStudio** (recommended for .Rmd notebooks)
- **Git** for cloning repository
- **Conda/Miniconda** (optional but recommended)

### Option 1: Using Conda (Recommended)

```bash
# Clone the repository
git clone https://github.com/yourusername/BCL2A1-LUAD-ICB-Biomarker.git
cd BCL2A1-LUAD-ICB-Biomarker

# Create conda environment
conda env create -f environment.yml
conda activate bcl2a1-analysis

# Install R packages
Rscript install_packages.R
```

### Option 2: Manual Installation

#### Step 1: Clone Repository

```bash
git clone https://github.com/yourusername/BCL2A1-LUAD-ICB-Biomarker.git
cd BCL2A1-LUAD-ICB-Biomarker
```

#### Step 2: Install R Packages

```bash
Rscript install_packages.R
```

Or manually in R:

```r
# CRAN packages
install.packages(c(
  "tidyverse", "data.table", "survival", "survminer",
  "caret", "glmnet", "pROC", "ggplot2", "pheatmap",
  "cowplot", "ggpubr", "ggrepel", "RColorBrewer",
  "reshape2", "scales", "gridExtra"
))

# Bioconductor packages
if (!require("BiocManager", quietly = TRUE))
    install.packages("BiocManager")

BiocManager::install(c(
  "DESeq2", "limma", "edgeR", "GSVA", "GSEABase",
  "fgsea", "clusterProfiler", "enrichplot",
  "org.Hs.eg.db", "GEOquery", "sva", "ComplexHeatmap"
))

# EPIC for deconvolution
devtools::install_github("GfellerLab/EPIC", build_vignettes=TRUE)
```

#### Step 3: Install Python Packages

```bash
pip install -r requirements.txt
```

### Verify Installation

```bash
# Test R installation
Rscript -e "library(tidyverse); library(DESeq2); library(EPIC)"

# Test Python installation
python -c "import scanpy; import pandas; import numpy"
```

---

## Data Availability

All datasets are publicly available from NCBI Gene Expression Omnibus (GEO):

### Bulk RNA-seq Cohorts

| Dataset | Description | Samples | Platform | Treatment | Reference |
|---------|-------------|---------|----------|-----------|-----------|
| **GSE161537** | Discovery cohort | 60 | Illumina NovaSeq | Anti-PD-(L)1 | Jung et al., 2021 |
| **GSE126044** | Validation cohort 1 | 7 | Illumina HiSeq | Anti-PD-1 | Cho et al., 2020 |
| **GSE135222** | Validation cohort 2 | 27 | Illumina HiSeq | Anti-PD-1 | Jung et al., 2019 |
| **GSE274975** | Validation cohort 3 | 28 | Illumina NovaSeq | Anti-PD-L1 | Recent, 2024 |

### Single-cell RNA-seq

| Dataset | Description | Patients | CD8+ Cells | Platform | Treatment |
|---------|-------------|----------|------------|----------|-----------|
| **GSE176021** | Neoadjuvant cohort | 11 | 130,584 | 10X Chromium | Neoadjuvant ICB + Chemo |

### Download Data

```bash
# Automated download (requires GEOquery)
cd data/
bash download_data.sh

# Manual download instructions in data/README.md
```

**Note:** Raw data files (>10GB) are not included in the repository. The download script fetches data directly from GEO.

---

## Quick Start

### Reproduce Main Figures

```r
# Open RStudio
# Run the figure reproduction notebook
rmarkdown::render("notebooks/03_reproduce_figures.Rmd")
```

This will generate all 8 main figures in `results/figures/`.

### Run Complete Pipeline

```bash
# Execute all modules sequentially
bash run_full_pipeline.sh
```

Or run individual modules:

```r
# Module 0: Preliminary Analysis
source("src/module_0_preliminary/01_differential_expression.R")
source("src/module_0_preliminary/02_survival_analysis.R")
# ... etc

# Module A: Discovery Modeling
source("src/module_A_discovery/01_data_preprocessing.R")
# ... etc
```

### Expected Runtime

| Module | Time | Memory | CPU |
|--------|------|--------|-----|
| Module 0 | ~30 min | 8 GB | 4 cores |
| Module A | ~2 hours | 16 GB | 8 cores |
| Module B | ~1 hour | 16 GB | 4 cores |
| Module C | ~3 hours | 32 GB | 8 cores |
| Module D | ~1 hour | 16 GB | 4 cores |
| Module E | ~45 min | 8 GB | 4 cores |
| **Total** | **~8 hours** | **32 GB** | **8 cores** |

---

## Module Descriptions

### Module 0: Preliminary Biomarker Characterization

**Purpose:** Identify and characterize BCL2A1 as a candidate biomarker

**Analyses:**
- Differential expression (responders vs non-responders)
- Survival analysis (Kaplan-Meier, Cox regression)
- Meta-analysis across 21 non-ICB LUAD cohorts
- Immune correlation with effector markers
- Pathway enrichment (GSEA, GO, KEGG)
- EPIC deconvolution (immune vs tumor compartments)

**Key Outputs:**
- Figure 1: BCL2A1 expression and survival
- Figure 2: Pathway enrichment and immune positioning
- Figure 3: Co-expression networks and deconvolution
- Supplementary Tables S1-S3

**Scripts:**
1. `01_differential_expression.R` - DE analysis with limma-voom
2. `02_survival_analysis.R` - KM curves and Cox models
3. `03_immune_correlation.R` - Correlation with immune markers
4. `04_pathway_enrichment.R` - GSEA and enrichment analysis
5. `05_meta_analysis.R` - Meta-analysis of non-ICB cohorts
6. `06_deconvolution_epic.R` - EPIC deconvolution

---

### Module A: Discovery Modeling

**Purpose:** Build and evaluate predictive models in discovery cohort

**Analyses:**
- Feature engineering (BCL2A1, CD274, HOT score)
- Logistic regression with L1 regularization
- Nested cross-validation (10 outer × 5 inner folds)
- Model comparison (with/without BCL2A1)
- Clinical augmentation (stage, treatment line)

**Key Outputs:**
- Figure 4: Discovery model performance
- Supplementary Tables S4-S8
- Model objects: `tri_marker_model.rds`

**Scripts:**
1. `01_data_preprocessing.R` - Normalization and filtering
2. `02_feature_engineering.R` - Compute HOT score
3. `03_model_training.R` - Train logistic regression models
4. `04_nested_cv.R` - Nested CV implementation
5. `05_model_evaluation.R` - Performance metrics
6. `06_comparative_analysis.R` - Model comparisons

---

### Module B: External Validation

**Purpose:** Validate models across independent cohorts

**Analyses:**
- Leave-one-cohort-out (LOCO) validation
- Batch correction (ComBat, quantile normalization)
- Model calibration (intercept-only recalibration)
- Decision curve analysis
- Clinical utility assessment

**Key Outputs:**
- Figure 5: External validation results
- Supplementary Tables S9-S11

**Scripts:**
1. `01_load_validation_cohorts.R` - Load GSE126044, GSE135222, GSE274975
2. `02_batch_correction.R` - Harmonize across cohorts
3. `03_loco_validation.R` - LOCO framework
4. `04_calibration.R` - Calibration analysis
5. `05_decision_curve_analysis.R` - DCA and net benefit

---

### Module C: Single-cell RNA-seq Analysis

**Purpose:** Characterize BCL2A1+ CD8+ T-cell states

**Analyses:**
- Quality control and filtering (130,584 cells)
- Normalization and scaling
- Dimensionality reduction (PCA, UMAP)
- Graph-based clustering (Leiden algorithm)
- CD8+ subtype annotation
- BCL2A1 expression profiling
- Differential abundance testing
- Subtype-specific marker extraction

**Key Outputs:**
- Figure 6: Single-cell CD8+ subtypes and BCL2A1 distribution
- Supplementary Table S12
- Signature matrices for deconvolution

**Scripts (Python):**
1. `01_data_loading.py` - Load GSE176021
2. `02_quality_control.py` - QC filtering
3. `03_normalization_scaling.py` - Preprocessing
4. `04_clustering.py` - Leiden clustering
5. `05_annotation.py` - Subtype annotation
6. `06_bcl2a1_analysis.py` - BCL2A1 expression analysis
7. `07_differential_abundance.py` - MPR vs non-MPR
8. `08_signature_extraction.py` - Extract marker genes

---

### Module D: Bulk Deconvolution & Robustness

**Purpose:** Validate single-cell findings in bulk transcriptomes

**Analyses:**
- Custom CD8+ reference matrix building
- EPIC deconvolution with subtype signatures
- Enrichment analysis (R vs NR)
- Ablation experiments (remove BCL2A1, correlated genes)
- Bootstrap resampling (1,000 iterations)
- Gaussian noise perturbation

**Key Outputs:**
- Figure 7: Deconvolution results and robustness
- Supplementary Tables S13-S17

**Scripts:**
1. `01_reference_building.R` - Build custom EPIC reference
2. `02_epic_deconvolution.R` - Run deconvolution
3. `03_enrichment_analysis.R` - Test subtype enrichment
4. `04_ablation_experiments.R` - Marker ablation
5. `05_robustness_testing.R` - Bootstrap and noise tests

---

### Module E: Cross-Platform Simulation

**Purpose:** Assess clinical translatability across platforms

**Analyses:**
- In silico simulation of 3 platforms:
  - HTG EdgeSeq
  - NanoString nCounter
  - RT-qPCR
- Platform-specific noise injection (CV, dropout, batch effects)
- Concordance analysis (Spearman ρ, Cohen's κ, Lin's CCC)
- Stress testing (increasing perturbations)
- Negative controls (single-gene models)

**Key Outputs:**
- Figure 8: Cross-platform concordance
- Supplementary Tables S18-S21
- Platform recommendations

**Scripts:**
1. `01_platform_simulation.R` - Simulate platform-specific data
2. `02_concordance_analysis.R` - Compute concordance metrics
3. `03_stress_testing.R` - Robustness under perturbations
4. `04_negative_controls.R` - Single-gene validation
5. `05_platform_comparison.R` - Platform rankings

---

## Reproducibility

### Computational Environment

The analysis was performed using:
- R version 4.3.3
- Python version 3.9.18
- RStudio 2023.12.1
- Ubuntu 22.04 LTS

All package versions are specified in:
- `requirements.txt` (Python)
- `install_packages.R` (R)
- `environment.yml` (Conda)

### Random Seeds

All analyses use fixed random seeds for reproducibility:
- R: `set.seed(42)`
- Python: `np.random.seed(42)`, `random.seed(42)`

### Session Info

Run `sessionInfo()` in R or check `docs/session_info.txt` for complete environment details.

---

## System Requirements

### Minimum Requirements
- **OS:** Linux, macOS, or Windows 10/11
- **RAM:** 16 GB
- **Storage:** 50 GB free space
- **CPU:** 4 cores @ 2.5 GHz

### Recommended Configuration
- **OS:** Linux (Ubuntu 20.04+) or macOS
- **RAM:** 32 GB (64 GB for Module C)
- **Storage:** 100 GB SSD
- **CPU:** 8+ cores @ 3.0 GHz
- **GPU:** Not required

---

## Troubleshooting

### Common Issues

#### 1. Package Installation Errors

**Problem:** Bioconductor packages fail to install

**Solution:**
```r
# Update Bioconductor
BiocManager::install(version = "3.18", update = TRUE, ask = FALSE)

# Install individual packages
BiocManager::install("EPIC", force = TRUE)
```

#### 2. Memory Errors

**Problem:** "Cannot allocate vector of size..."

**Solution:**
```r
# Increase memory limit (Windows)
memory.limit(size = 32000)

# For Linux/Mac, run R with:
R --max-mem-size=32G
```

#### 3. GEO Download Timeouts

**Problem:** Timeout errors when downloading data

**Solution:**
```r
options(timeout = 600)
library(GEOquery)
gse <- getGEO("GSE161537", GSEMatrix = TRUE, AnnotGPL = FALSE)
```

#### 4. Python scanpy Issues

**Problem:** `ImportError: cannot import name 'read' from 'scanpy'`

**Solution:**
```bash
pip install --upgrade scanpy anndata
```

For more solutions, see `docs/troubleshooting.md`.

---

## Contributing

We welcome contributions! To contribute:

1. Fork the repository
2. Create a feature branch (`git checkout -b feature/amazing-feature`)
3. Commit your changes (`git commit -m 'Add amazing feature'`)
4. Push to the branch (`git push origin feature/amazing-feature`)
5. Open a Pull Request

Please ensure:
- Code follows existing style conventions
- New functions include documentation
- Tests pass (if applicable)
- Update README.md if needed

---

## Contact

**Corresponding Author:**  
**Prof. Chiou-Feng Lin**  
Department of Microbiology and Immunology  
School of Medicine, College of Medicine  
Taipei Medical University  
No. 250, Wuxing Street, Taipei 11031, Taiwan  
📧 Email: cflin2014@tmu.edu.tw  
🌐 Lab Website: [Link]

**First Author:**  
**Hoang Minh Quan Pham**  
International Ph.D. Program in Medicine  
Taipei Medical University  
📧 Email: [your.email@tmu.edu.tw]

**Issues & Questions:**  
Please report bugs or request features via [GitHub Issues](https://github.com/yourusername/BCL2A1-LUAD-ICB-Biomarker/issues)

---

## License

This project is licensed under the **MIT License** - see the [LICENSE](LICENSE) file for details.

```
MIT License

Copyright (c) 2025 Hoang Minh Quan Pham, Chiou-Feng Lin

Permission is hereby granted, free of charge, to any person obtaining a copy
of this software and associated documentation files (the "Software"), to deal
in the Software without restriction...
```

---

## Acknowledgments

### Funding

This work was funded by grants from the **National Science and Technology Council (NSTC)**, Taiwan:
- NSTC113-2314-B-038-137
- NSTC114-2314-B-038-013

### Data Contributors

We thank the authors of the following datasets for making their data publicly available:
- Jung et al. (GSE161537, GSE135222)
- Cho et al. (GSE126044)
- GSE274975 contributors
- GSE176021 contributors

### Software & Tools

This analysis builds upon excellent open-source tools:
- [DESeq2](https://bioconductor.org/packages/DESeq2/) (Love et al.)
- [Scanpy](https://scanpy.readthedocs.io/) (Wolf et al.)
- [EPIC](https://github.com/GfellerLab/EPIC) (Racle et al.)
- [clusterProfiler](https://bioconductor.org/packages/clusterProfiler/) (Yu et al.)
- And many others listed in requirements

---

## Version History

- **v1.0.0** (2025-01-XX): Initial release with manuscript submission to CSBJ
- **v1.0.1** (TBD): Bug fixes and minor improvements

---

## Citation Statistics

**BibTeX:**
```bibtex
@article{pham2025bcl2a1,
  title={Integrative Computational Analysis Identifies BCL2A1 as a CD8+ T-cell Survival Marker Associated with Immunotherapy Response in Lung Adenocarcinoma},
  author={Pham, Hoang Minh Quan and Feng, Po-Hao and Chen, Chia-Ling and Lee, Kang-Yun and Lin, Chiou-Feng},
  journal={Computational and Structural Biotechnology Journal},
  year={2025},
  doi={10.XXXX/XXXXX},
  url={https://github.com/yourusername/BCL2A1-LUAD-ICB-Biomarker}
}
```

**APA:**
Pham, H. M. Q., Feng, P.-H., Chen, C.-L., Lee, K.-Y., & Lin, C.-F. (2025). Integrative computational analysis identifies BCL2A1 as a CD8⁺ T-cell survival marker associated with immunotherapy response in lung adenocarcinoma. *Computational and Structural Biotechnology Journal*. https://doi.org/10.XXXX/XXXXX

---

**Last Updated:** January 2025  
**Repository Status:** ✅ Active Development

For the latest updates, visit: https://github.com/yourusername/BCL2A1-LUAD-ICB-Biomarker

---

*Made with ❤️ for reproducible science*