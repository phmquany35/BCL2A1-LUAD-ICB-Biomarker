#!/usr/bin/env Rscript

#' @title Plotting Utilities
#' @description Common plotting functions for visualization
#' @author Hoang Minh Quan Pham

suppressPackageStartupMessages({
  library(ggplot2)
  library(ggrepel)
  library(ggpubr)
  library(pheatmap)
  library(RColorBrewer)
  library(cowplot)
  library(scales)
})

# Set default theme
theme_set(theme_bw() + 
          theme(panel.grid.minor = element_blank(),
                axis.text = element_text(size = 10),
                axis.title = element_text(size = 12),
                plot.title = element_text(size = 14, face = "bold")))

#' Volcano plot for differential expression
#' 
#' @param results DE results data frame
#' @param logfc_col Column name for log fold change
#' @param pval_col Column name for p-value
#' @param fdr_col Column name for adjusted p-value
#' @param fdr_threshold FDR threshold (default: 0.05)
#' @param logfc_threshold Log FC threshold (default: 1)
#' @param label_top Number of top genes to label (default: 10)
#' @return ggplot object
plot_volcano <- function(results, 
                        logfc_col = "logFC",
                        pval_col = "P.Value",
                        fdr_col = "adj.P.Val",
                        fdr_threshold = 0.05,
                        logfc_threshold = 1,
                        label_top = 10) {
  
  # Prepare data
  plot_data <- data.frame(
    gene = rownames(results),
    logFC = results[[logfc_col]],
    pvalue = results[[pval_col]],
    fdr = results[[fdr_col]]
  )
  
  # Add significance category
  plot_data$sig <- "NS"
  plot_data$sig[plot_data$fdr < fdr_threshold & plot_data$logFC > logfc_threshold] <- "Up"
  plot_data$sig[plot_data$fdr < fdr_threshold & plot_data$logFC < -logfc_threshold] <- "Down"
  plot_data$sig <- factor(plot_data$sig, levels = c("Down", "NS", "Up"))
  
  # Select genes to label
  sig_genes <- plot_data[plot_data$sig != "NS", ]
  sig_genes <- sig_genes[order(sig_genes$pvalue), ]
  genes_to_label <- head(sig_genes, label_top)
  
  # Plot
  p <- ggplot(plot_data, aes(x = logFC, y = -log10(pvalue), color = sig)) +
    geom_point(alpha = 0.6, size = 1.5) +
    geom_text_repel(
      data = genes_to_label,
      aes(label = gene),
      size = 3,
      max.overlaps = 20
    ) +
    scale_color_manual(
      values = c("Down" = "#2166AC", "NS" = "grey70", "Up" = "#B2182B"),
      labels = c(
        sprintf("Down (n=%d)", sum(plot_data$sig == "Down")),
        sprintf("NS (n=%d)", sum(plot_data$sig == "NS")),
        sprintf("Up (n=%d)", sum(plot_data$sig == "Up"))
      )
    ) +
    geom_hline(yintercept = -log10(0.05), linetype = "dashed", color = "grey30") +
    geom_vline(xintercept = c(-logfc_threshold, logfc_threshold), 
               linetype = "dashed", color = "grey30") +
    labs(
      title = "Volcano Plot: Responders vs Non-Responders",
      x = expression(Log[2]~Fold~Change),
      y = expression(-Log[10]~P-value),
      color = "Significance"
    ) +
    theme(legend.position = "top")
  
  return(p)
}


#' Kaplan-Meier survival plot
#'
#' @param survfit_object survfit object from survival package
#' @param title Plot title
#' @param xlab X-axis label
#' @param ylab Y-axis label
#' @return ggsurvplot object
plot_km_curve <- function(survfit_object,
                         title = "Kaplan-Meier Survival Curve",
                         xlab = "Time (months)",
                         ylab = "Overall Survival Probability") {
  
  require(survminer)
  
  p <- ggsurvplot(
    survfit_object,
    data = survfit_object$call$data,
    pval = TRUE,
    conf.int = TRUE,
    risk.table = TRUE,
    risk.table.height = 0.25,
    ggtheme = theme_bw(),
    palette = c("#E7B800", "#2E9FDF"),
    title = title,
    xlab = xlab,
    ylab = ylab,
    legend.title = "",
    legend.labs = c("Low", "High"),
    break.x.by = 6
  )
  
  return(p)
}


#' Heatmap with clustering
#'
#' @param mat Matrix for heatmap
#' @param annotation_col Column annotations
#' @param scale Scale rows, columns, or none
#' @param cluster_rows Cluster rows
#' @param cluster_cols Cluster columns
#' @param show_rownames Show row names
#' @param show_colnames Show column names
#' @param title Plot title
#' @return pheatmap object
plot_heatmap <- function(mat,
                        annotation_col = NULL,
                        scale = "row",
                        cluster_rows = TRUE,
                        cluster_cols = TRUE,
                        show_rownames = FALSE,
                        show_colnames = FALSE,
                        title = "") {
  
  # Color palette
  col_palette <- colorRampPalette(rev(brewer.pal(11, "RdBu")))(100)
  
  # Plot
  p <- pheatmap(
    mat,
    color = col_palette,
    scale = scale,
    cluster_rows = cluster_rows,
    cluster_cols = cluster_cols,
    show_rownames = show_rownames,
    show_colnames = show_colnames,
    annotation_col = annotation_col,
    main = title,
    fontsize = 10,
    fontsize_row = 8,
    fontsize_col = 8
  )
  
  return(p)
}


#' Boxplot comparing groups
#'
#' @param data Data frame
#' @param x_var X-axis variable (grouping)
#' @param y_var Y-axis variable (values)
#' @param fill_var Fill variable (optional)
#' @param title Plot title
#' @param ylab Y-axis label
#' @return ggplot object
plot_boxplot <- function(data, x_var, y_var, fill_var = NULL, 
                        title = "", ylab = "Expression") {
  
  p <- ggplot(data, aes_string(x = x_var, y = y_var)) +
    {if (!is.null(fill_var)) geom_boxplot(aes_string(fill = fill_var))
     else geom_boxplot(fill = "steelblue")} +
    geom_jitter(width = 0.2, alpha = 0.5, size = 2) +
    stat_compare_means() +
    labs(title = title, y = ylab) +
    theme(axis.text.x = element_text(angle = 45, hjust = 1))
  
  return(p)
}


#' ROC curve plot
#'
#' @param roc_obj ROC object from pROC package
#' @param title Plot title
#' @return ggplot object
plot_roc_curve <- function(roc_obj, title = "ROC Curve") {
  
  roc_data <- data.frame(
    fpr = 1 - roc_obj$specificities,
    tpr = roc_obj$sensitivities
  )
  
  auc_value <- round(auc(roc_obj), 3)
  
  p <- ggplot(roc_data, aes(x = fpr, y = tpr)) +
    geom_line(color = "#2E9FDF", size = 1.2) +
    geom_abline(linetype = "dashed", color = "grey50") +
    annotate("text", x = 0.7, y = 0.3, 
             label = sprintf("AUC = %.3f", auc_value),
             size = 5, fontface = "bold") +
    labs(
      title = title,
      x = "False Positive Rate (1 - Specificity)",
      y = "True Positive Rate (Sensitivity)"
    ) +
    coord_fixed() +
    theme_minimal()
  
  return(p)
}


#' Save plot to file
#'
#' @param plot ggplot or pheatmap object
#' @param filename Output filename
#' @param width Width in inches
#' @param height Height in inches
#' @param dpi Resolution
save_plot_file <- function(plot, filename, width = 8, height = 6, dpi = 300) {
  
  # Create directory if doesn't exist
  output_dir <- dirname(filename)
  if (!dir.exists(output_dir)) {
    dir.create(output_dir, recursive = TRUE)
  }
  
  # Save
  ggsave(filename, plot, width = width, height = height, dpi = dpi)
  cat(sprintf("Plot saved to %s\n", filename))
}

cat("✓ plotting_functions.R loaded\n")
