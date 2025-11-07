#!/usr/bin/env Rscript
# Differential Expression Analysis
# Author: Hoang Minh Quan Pham

suppressPackageStartupMessages({
  library(tidyverse)
  library(limma)
  library(edgeR)
})

source("src/utils/data_loader.R")
source("src/utils/plotting_functions.R")

main <- function() {
  
  cat("="*50, "\n")
  cat("MODULE 0.1: Differential Expression Analysis\n")
  cat("="*50, "\n\n")
  
  # Load discovery cohort
  data <- load_discovery_cohort()
  expr <- data$expression
  meta <- data$metadata
  
  # Design matrix
  design <- model.matrix(~0 + Response, data = meta)
  colnames(design) <- gsub("Response", "", colnames(design))
  
  # Voom transformation
  v <- voom(expr, design, plot = FALSE)
  
  # Fit linear model
  fit <- lmFit(v, design)
  
  # Contrast
  contrast.matrix <- makeContrasts(
    RvsNR = Responder - `Non-Responder`,
    levels = design
  )
  
  fit2 <- contrasts.fit(fit, contrast.matrix)
  fit2 <- eBayes(fit2)
  
  # Results
  results <- topTable(fit2, coef = "RvsNR", number = Inf)
  
  # Save
  write.csv(results, "results/tables/differential_expression_all.csv")
  
  # Significant genes
  sig_genes <- results %>%
    filter(adj.P.Val < 0.05) %>%
    arrange(desc(abs(logFC)))
  
  write.csv(sig_genes, "results/tables/differential_expression_significant.csv")
  
  cat(sprintf("Found %d significant genes\n", nrow(sig_genes)))
  
  # Volcano plot
  p <- plot_volcano(results)
  save_plot_file(p, "results/figures/Figure_1A_volcano.pdf")
  
  cat("\nModule 0.1 complete!\n")
  
  return(list(results = results, significant = sig_genes))
}

if (!interactive()) {
  results <- main()
}
