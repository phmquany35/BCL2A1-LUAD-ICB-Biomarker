#!/usr/bin/env Rscript
# Model Training - Tri-Marker Model
# Author: Hoang Minh Quan Pham

suppressPackageStartupMessages({
  library(glmnet)
  library(caret)
  library(pROC)
  library(tidyverse)
})

source("src/utils/data_loader.R")

main <- function() {
  
  cat("MODULE A.3: Model Training\n\n")
  
  # Load data
  data <- load_discovery_cohort()
  expr <- data$expression
  meta <- data$metadata
  
  # Features: BCL2A1, CD274, HOT score
  features <- data.frame(
    BCL2A1 = expr["BCL2A1", ],
    CD274 = expr["CD274", ],
    HOT_score = colMeans(expr[hot_genes, ])  # Define hot_genes from signature
  )
  
  # Outcome
  y <- as.numeric(meta$Response == "Responder")
  
  # Train logistic regression with L1 penalty
  set.seed(42)
  cv_fit <- cv.glmnet(as.matrix(features), y, family = "binomial", alpha = 1)
  
  # Predictions
  pred_prob <- predict(cv_fit, newx = as.matrix(features), 
                      s = "lambda.min", type = "response")
  
  # ROC analysis
  roc_obj <- roc(y, as.vector(pred_prob))
  cat(sprintf("AUC = %.3f\n", auc(roc_obj)))
  
  # Save model
  saveRDS(cv_fit, "results/models/tri_marker_model.rds")
  
  # Plot ROC
  p <- plot_roc_curve(roc_obj, "Tri-Marker Model Performance")
  save_plot_file(p, "results/figures/Figure_4D_roc.pdf")
  
  cat("\nModule A.3 complete!\n")
}

if (!interactive()) {
  main()
}
