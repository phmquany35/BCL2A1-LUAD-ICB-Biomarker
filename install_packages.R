#!/usr/bin/env Rscript

#' @title Install R Packages for BCL2A1 Analysis
#' @description Automated installation of all required R packages
#' @author Hoang Minh Quan Pham
#' @date 2025-01-XX

cat("====================================\n")
cat("BCL2A1 Analysis Package Installer\n")
cat("====================================\n\n")

# Set CRAN mirror
options(repos = c(CRAN = "https://cloud.r-project.org"))

# Function to install packages if not already installed
install_if_missing <- function(packages, type = "CRAN") {
  for (pkg in packages) {
    if (!require(pkg, character.only = TRUE, quietly = TRUE)) {
      cat(sprintf("Installing %s from %s...\n", pkg, type))
      if (type == "CRAN") {
        install.packages(pkg, dependencies = TRUE)
      } else if (type == "Bioconductor") {
        BiocManager::install(pkg, update = FALSE, ask = FALSE)
      } else if (type == "GitHub") {
        devtools::install_github(pkg, upgrade = "never")
      }
    } else {
      cat(sprintf("✓ %s already installed\n", pkg))
    }
  }
}

# ==========================================
# 1. Install BiocManager first
# ==========================================
cat("\n[1/5] Installing BiocManager...\n")
if (!require("BiocManager", quietly = TRUE)) {
  install.packages("BiocManager")
}
library(BiocManager)

# ==========================================
# 2. Install CRAN packages
# ==========================================
cat("\n[2/5] Installing CRAN packages...\n")

cran_packages <- c(
  # Data manipulation
  "tidyverse",
  "data.table",
  "reshape2",
  "dplyr",
  "tidyr",
  "tibble",
  
  # Statistics
  "survival",
  "survminer",
  "caret",
  "glmnet",
  "pROC",
  "boot",
  "meta",
  "metafor",
  
  # Visualization
  "ggplot2",
  "ggrepel",
  "ggpubr",
  "pheatmap",
  "cowplot",
  "RColorBrewer",
  "scales",
  "gridExtra",
  "viridis",
  "ggsci",
  
  # Utilities
  "devtools",
  "remotes",
  "here",
  "readr",
  "writexl",
  "openxlsx"
)

install_if_missing(cran_packages, type = "CRAN")

# ==========================================
# 3. Install Bioconductor packages
# ==========================================
cat("\n[3/5] Installing Bioconductor packages...\n")

bioc_packages <- c(
  # RNA-seq analysis
  "DESeq2",
  "limma",
  "edgeR",
  
  # Gene set analysis
  "GSVA",
  "GSEABase",
  "fgsea",
  "clusterProfiler",
  "enrichplot",
  "msigdbr",
  
  # Annotation
  "org.Hs.eg.db",
  "AnnotationDbi",
  
  # Data retrieval
  "GEOquery",
  "Biobase",
  
  # Batch correction
  "sva",
  
  # Visualization
  "ComplexHeatmap",
  "circlize"
)

install_if_missing(bioc_packages, type = "Bioconductor")

# ==========================================
# 4. Install GitHub packages
# ==========================================
cat("\n[4/5] Installing GitHub packages...\n")

if (!require("devtools", quietly = TRUE)) {
  install.packages("devtools")
}
library(devtools)

# EPIC for deconvolution
if (!require("EPIC", quietly = TRUE)) {
  cat("Installing EPIC from GitHub...\n")
  devtools::install_github("GfellerLab/EPIC", build_vignettes = FALSE, upgrade = "never")
} else {
  cat("✓ EPIC already installed\n")
}

# ==========================================
# 5. Verify installation
# ==========================================
cat("\n[5/5] Verifying installation...\n")

essential_packages <- c(
  "tidyverse", "DESeq2", "limma", "survival", 
  "caret", "glmnet", "pROC", "ggplot2",
  "clusterProfiler", "EPIC", "GEOquery", "sva"
)

all_installed <- TRUE
for (pkg in essential_packages) {
  if (require(pkg, character.only = TRUE, quietly = TRUE)) {
    cat(sprintf("✓ %s loaded successfully\n", pkg))
  } else {
    cat(sprintf("✗ %s failed to load\n", pkg))
    all_installed <- FALSE
  }
}

# ==========================================
# Summary
# ==========================================
cat("\n====================================\n")
if (all_installed) {
  cat("✓ All packages installed successfully!\n")
  cat("====================================\n")
  cat("\nYou can now run the analysis pipeline.\n")
  cat("To get started, try:\n")
  cat("  source('src/module_0_preliminary/01_differential_expression.R')\n\n")
} else {
  cat("✗ Some packages failed to install.\n")
  cat("====================================\n")
  cat("\nPlease check the error messages above and try:\n")
  cat("  1. Update R to version >= 4.3.0\n")
  cat("  2. Update Bioconductor: BiocManager::install(version='3.18')\n")
  cat("  3. Install failed packages manually\n\n")
}

# Print session info
cat("\nR Session Information:\n")
cat("----------------------\n")
print(sessionInfo())
