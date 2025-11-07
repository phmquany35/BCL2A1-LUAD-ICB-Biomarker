#!/usr/bin/env Rscript
# Survival Analysis
# Author: Hoang Minh Quan Pham

suppressPackageStartupMessages({
  library(survival)
  library(survminer)
  library(tidyverse)
})

source("src/utils/data_loader.R")
source("src/utils/plotting_functions.R")

main <- function() {
  
  cat("MODULE 0.2: Survival Analysis\n\n")
  
  # Load data
  data <- load_discovery_cohort()
  expr <- data$expression
  meta <- data$metadata
  
  # BCL2A1 expression
  bcl2a1_expr <- expr["BCL2A1", ]
  meta$BCL2A1 <- bcl2a1_expr[rownames(meta)]
  meta$BCL2A1_group <- ifelse(meta$BCL2A1 > median(meta$BCL2A1), "High", "Low")
  
  # Survival analysis
  surv_obj <- Surv(time = meta$OS_months, event = meta$OS_event)
  fit <- survfit(surv_obj ~ BCL2A1_group, data = meta)
  
  # Plot
  p <- plot_km_curve(fit, title = "BCL2A1 Expression and Overall Survival")
  pdf("results/figures/Figure_1C_survival.pdf", width = 8, height = 6)
  print(p)
  dev.off()
  
  # Cox regression
  cox_model <- coxph(surv_obj ~ BCL2A1 + Response, data = meta)
  cat("\nCox Regression Results:\n")
  print(summary(cox_model))
  
  cat("\nModule 0.2 complete!\n")
}

if (!interactive()) {
  main()
}
