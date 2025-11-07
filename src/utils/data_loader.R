#!/usr/bin/env Rscript

#' @title Data Loading Utilities
#' @description Functions for loading and preprocessing GEO datasets
#' @author Hoang Minh Quan Pham
#' @date 2025-01-XX

# Load required libraries
suppressPackageStartupMessages({
  library(GEOquery)
  library(Biobase)
  library(tidyverse)
})

#' Load GEO dataset with expression matrix and metadata
#'
#' @param geo_id GEO accession ID (e.g., "GSE161537")
#' @param data_dir Directory to save/load data (default: "data/raw/")
#' @param force_download Force re-download even if file exists
#' @return List containing expression matrix and metadata
#' @export
#'
#' @examples
#' data <- load_geo_data("GSE161537")
#' expr <- data$expression
#' meta <- data$metadata
load_geo_data <- function(geo_id, 
                         data_dir = "data/raw/", 
                         force_download = FALSE) {
  
  cat(sprintf("Loading GEO dataset: %s\n", geo_id))
  
  # Create data directory if doesn't exist
  if (!dir.exists(data_dir)) {
    dir.create(data_dir, recursive = TRUE)
    cat(sprintf("Created directory: %s\n", data_dir))
  }
  
  # Set options for download
  options(timeout = 600)  # 10 minutes timeout
  
  # Check if already downloaded
  rds_file <- file.path(data_dir, paste0(geo_id, ".rds"))
  
  if (file.exists(rds_file) && !force_download) {
    cat(sprintf("Loading cached data from %s\n", rds_file))
    data <- readRDS(rds_file)
    return(data)
  }
  
  # Download from GEO
  cat("Downloading from GEO (this may take several minutes)...\n")
  
  tryCatch({
    # Get GEO data
    gse <- getGEO(geo_id, GSEMatrix = TRUE, AnnotGPL = FALSE, destdir = data_dir)
    
    # Extract expression set
    if (length(gse) > 1) {
      warning("Multiple platforms found. Using the first one.")
      eset <- gse[[1]]
    } else {
      eset <- gse[[1]]
    }
    
    # Extract expression matrix
    expr_matrix <- exprs(eset)
    cat(sprintf("Expression matrix: %d genes × %d samples\n", 
                nrow(expr_matrix), ncol(expr_matrix)))
    
    # Extract metadata
    metadata <- pData(eset)
    cat(sprintf("Metadata: %d samples × %d variables\n", 
                nrow(metadata), ncol(metadata)))
    
    # Extract feature data (gene annotations)
    feature_data <- fData(eset)
    
    # Prepare output
    output <- list(
      expression = expr_matrix,
      metadata = metadata,
      features = feature_data,
      eset = eset
    )
    
    # Save to RDS for faster loading next time
    saveRDS(output, rds_file)
    cat(sprintf("Saved data to %s\n", rds_file))
    
    return(output)
    
  }, error = function(e) {
    cat(sprintf("Error downloading %s: %s\n", geo_id, e$message))
    cat("\nTroubleshooting tips:\n")
    cat("1. Check internet connection\n")
    cat("2. Verify GEO ID is correct\n")
    cat("3. Try manual download from https://www.ncbi.nlm.nih.gov/geo/\n")
    cat("4. Increase timeout: options(timeout=900)\n")
    return(NULL)
  })
}


#' Load multiple GEO datasets
#'
#' @param geo_ids Vector of GEO accession IDs
#' @param data_dir Directory to save/load data
#' @return Named list of datasets
#' @export
load_multiple_geo_data <- function(geo_ids, data_dir = "data/raw/") {
  
  cat(sprintf("Loading %d GEO datasets...\n", length(geo_ids)))
  
  datasets <- list()
  for (geo_id in geo_ids) {
    cat(sprintf("\n[%d/%d] ", which(geo_ids == geo_id), length(geo_ids)))
    datasets[[geo_id]] <- load_geo_data(geo_id, data_dir)
  }
  
  cat("\nAll datasets loaded successfully!\n")
  return(datasets)
}


#' Preprocess expression matrix
#'
#' @param expr_matrix Expression matrix (genes × samples)
#' @param log2_transform Apply log2 transformation
#' @param filter_low_expr Filter genes with low expression
#' @param min_expr_threshold Minimum expression threshold
#' @param min_sample_prop Minimum proportion of samples with expression
#' @return Preprocessed expression matrix
#' @export
preprocess_expression <- function(expr_matrix,
                                  log2_transform = TRUE,
                                  filter_low_expr = TRUE,
                                  min_expr_threshold = 1,
                                  min_sample_prop = 0.1) {
  
  cat("Preprocessing expression matrix...\n")
  cat(sprintf("Input: %d genes × %d samples\n", nrow(expr_matrix), ncol(expr_matrix)))
  
  # Filter low expression genes
  if (filter_low_expr) {
    min_samples <- ceiling(ncol(expr_matrix) * min_sample_prop)
    keep <- rowSums(expr_matrix > min_expr_threshold) >= min_samples
    expr_matrix <- expr_matrix[keep, ]
    cat(sprintf("After filtering: %d genes retained\n", nrow(expr_matrix)))
  }
  
  # Log2 transformation
  if (log2_transform) {
    # Check if already log-transformed
    max_val <- max(expr_matrix, na.rm = TRUE)
    if (max_val > 50) {
      cat("Applying log2 transformation...\n")
      expr_matrix <- log2(expr_matrix + 1)
    } else {
      cat("Data appears already log-transformed (skipping)\n")
    }
  }
  
  # Remove rows with NA or Inf
  expr_matrix <- expr_matrix[complete.cases(expr_matrix), ]
  expr_matrix <- expr_matrix[!apply(expr_matrix, 1, function(x) any(is.infinite(x))), ]
  
  cat(sprintf("Final: %d genes × %d samples\n", nrow(expr_matrix), ncol(expr_matrix)))
  
  return(expr_matrix)
}


#' Parse clinical metadata
#'
#' @param metadata Raw metadata from GEO
#' @param response_column Column name for response status
#' @param response_positive Positive response indicators
#' @param response_negative Negative response indicators
#' @return Cleaned metadata with binary response
#' @export
parse_clinical_metadata <- function(metadata,
                                   response_column = "response:ch1",
                                   response_positive = c("R", "PR", "CR", "MPR"),
                                   response_negative = c("NR", "PD", "SD", "non-MPR")) {
  
  cat("Parsing clinical metadata...\n")
  
  # Create clean copy
  meta <- as.data.frame(metadata)
  
  # Parse response status
  if (response_column %in% colnames(meta)) {
    response_raw <- meta[[response_column]]
    
    # Binarize response
    meta$Response <- NA
    meta$Response[response_raw %in% response_positive] <- "Responder"
    meta$Response[response_raw %in% response_negative] <- "Non-Responder"
    meta$Response <- factor(meta$Response, levels = c("Non-Responder", "Responder"))
    
    cat(sprintf("Response status parsed: %d responders, %d non-responders\n",
                sum(meta$Response == "Responder", na.rm = TRUE),
                sum(meta$Response == "Non-Responder", na.rm = TRUE)))
  } else {
    warning(sprintf("Response column '%s' not found in metadata", response_column))
  }
  
  return(meta)
}


#' Load preprocessed discovery cohort (GSE161537)
#'
#' @param data_dir Data directory
#' @return List with expression and metadata
#' @export
load_discovery_cohort <- function(data_dir = "data/raw/") {
  
  cat("Loading discovery cohort (GSE161537)...\n")
  
  # Load data
  data <- load_geo_data("GSE161537", data_dir)
  
  if (is.null(data)) {
    stop("Failed to load discovery cohort")
  }
  
  # Preprocess expression
  expr <- preprocess_expression(data$expression)
  
  # Parse metadata
  meta <- parse_clinical_metadata(data$metadata)
  
  # Match samples
  common_samples <- intersect(colnames(expr), rownames(meta))
  expr <- expr[, common_samples]
  meta <- meta[common_samples, ]
  
  cat(sprintf("Final dataset: %d genes × %d samples\n", 
              nrow(expr), ncol(expr)))
  
  return(list(
    expression = expr,
    metadata = meta,
    features = data$features
  ))
}


#' Load validation cohorts
#'
#' @param data_dir Data directory
#' @return Named list of validation cohorts
#' @export
load_validation_cohorts <- function(data_dir = "data/raw/") {
  
  validation_ids <- c("GSE126044", "GSE135222", "GSE274975")
  
  cat("Loading validation cohorts...\n")
  datasets <- load_multiple_geo_data(validation_ids, data_dir)
  
  # Preprocess each dataset
  for (geo_id in names(datasets)) {
    cat(sprintf("\nPreprocessing %s...\n", geo_id))
    
    if (!is.null(datasets[[geo_id]])) {
      datasets[[geo_id]]$expression <- preprocess_expression(datasets[[geo_id]]$expression)
      datasets[[geo_id]]$metadata <- parse_clinical_metadata(datasets[[geo_id]]$metadata)
    }
  }
  
  return(datasets)
}


#' Save processed data
#'
#' @param data Data object to save
#' @param filename Output filename
#' @param output_dir Output directory
#' @export
save_processed_data <- function(data, filename, output_dir = "results/") {
  
  if (!dir.exists(output_dir)) {
    dir.create(output_dir, recursive = TRUE)
  }
  
  filepath <- file.path(output_dir, filename)
  saveRDS(data, filepath)
  cat(sprintf("Data saved to %s\n", filepath))
}


#' Load processed data
#'
#' @param filename Input filename
#' @param input_dir Input directory
#' @return Loaded data object
#' @export
load_processed_data <- function(filename, input_dir = "results/") {
  
  filepath <- file.path(input_dir, filename)
  
  if (!file.exists(filepath)) {
    stop(sprintf("File not found: %s", filepath))
  }
  
  data <- readRDS(filepath)
  cat(sprintf("Data loaded from %s\n", filepath))
  
  return(data)
}


# Print message when loaded
cat("✓ data_loader.R loaded successfully\n")
cat("Available functions:\n")
cat("  - load_geo_data()\n")
cat("  - load_multiple_geo_data()\n")
cat("  - preprocess_expression()\n")
cat("  - parse_clinical_metadata()\n")
cat("  - load_discovery_cohort()\n")
cat("  - load_validation_cohorts()\n")
cat("  - save_processed_data()\n")
cat("  - load_processed_data()\n")
