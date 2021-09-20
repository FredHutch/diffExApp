## COUNTS TABLE VISUALIZATIONS ####################################################################
countsToPca <- function(de_out,
                        sample_matrix,
                        de_package) {
  # get counts from de_out
  if (de_package == "DESeq2") {
    counts <- counts(de_out, normalized = TRUE)
  } else if (de_package == "edgeR") {
  counts <- de_out$counts
  }
  pca <- prcomp(t(counts))
  pca_df <- data.frame(pca$x)
  pca_df$condition <- as.character(sample_matrix$condition)
  
  ggplot(pca_df, aes(x = PC1, y = PC2, color = condition)) +
    geom_point() +
    ggrepel::geom_label_repel(aes(label = rownames(pca_df)), show.legend = FALSE) +
    theme_classic()
}

## RESULTS VISUALIZATIONS #########################################################################
# plot de output in MA plot
resultsToMa <- function(de_out,
                        de_package,
                        pvalue_threshold,
                        logfc_threshold,
                        fdr) {
  # pull results from de output object
  if (de_package == "DESeq2") {
    de_res <- data.frame(results(de_out))
  } else if (de_package == "edgeR") {
    de_res <- data.frame(topTags(de_out, n = Inf))
  }
  # handle differing results table naming schemes
  if (de_package == "DESeq2") {
    mean_counts_col <- "baseMean"
    logfc_col <- "log2FoldChange"
    pvalue_col <- ifelse(fdr, "padj", "pvalue")

  } else if (de_package == "edgeR") {
    mean_counts_col <- "logCPM"
    logfc_col <- "logFC"
    pvalue_col <- ifelse(fdr, "FDR", "PValue")
  }
  
  # create boolean vector of significance
  de_idx <- abs(de_res[[logfc_col]]) >= logfc_threshold & de_res[[pvalue_col]] <= pvalue_threshold
  # NA to FALSE
  de_idx[is.na(de_idx)] = FALSE
  # Reset levels so TRUE is first
  de_idx <- factor(de_idx, levels = c(TRUE, FALSE))
  
  # plot
  ma <- ggplot(de_res, aes(.data[[mean_counts_col]], .data[[logfc_col]], colour = de_idx)) +
    geom_point(size=2, alpha = .6) + 
    scale_y_continuous(limits=c(-3, 3), oob = squish)
  
  # if package is DESeq log10 the base mean axis
  if (de_package == "DESeq2") {
    ma <- ma +
      scale_x_log10()
  }
  
  ma +
    geom_hline(yintercept = 0, colour="grey", size=1, linetype = "dashed") + 
    labs(x="Mean of Normalized Counts", y="Log Fold Change") + 
    theme_classic()
}

# plot de results as volcano plot
resultsToVolcano <- function(de_out,
                             de_package,
                             pvalue_threshold,
                             logfc_threshold,
                             fdr) {
  # pull results from de results object
  if (de_package == "DESeq2") {
    de_res <- data.frame(results(de_out))
  } else if (de_package == "edgeR") {
    de_res <- data.frame(topTags(de_out, n = Inf))
  }
  
  # handle differing results table naming schemes
  if (de_package == "DESeq2") {
    logfc_col <- "log2FoldChange"
    pvalue_col <- ifelse(fdr, "padj", "pvalue")
    
  } else if (de_package == "edgeR") {
    logfc_col <- "logFC"
    pvalue_col <- ifelse(fdr, "FDR", "PValue")
  }
  
  # convert pval to -log10(pval)
  de_res <- mutate(de_res,
                 log_pval = -log10(de_res[[pvalue_col]]))
  
  de_idx <- abs(de_res[[logfc_col]]) >= logfc_threshold & de_res[[pvalue_col]] <= pvalue_threshold
  # NA to FALSE
  de_idx[is.na(de_idx)] = FALSE
  # Reset levels so TRUE is first
  de_idx <- factor(de_idx, levels = c(TRUE, FALSE))

  # build base of plot
  volcano <- ggplot(de_res, aes(x = .data[[logfc_col]], y = log_pval)) +
    geom_point(size = 2, alpha = .6, aes(color = de_idx)) +
    geom_hline(yintercept = -log10(pvalue_threshold), linetype = "dashed", col = "grey", size = 1) +
    geom_vline(xintercept = c(logfc_threshold, -logfc_threshold), linetype = "dashed", col = "grey", size = 1) +
    scale_y_continuous(limits=c(0, 10), oob = squish) +
    labs(x = "Log Fold Change", y = "-log10(P value)", color = "Is DE") +
    theme_classic(base_size = 12)
  
  volcano

}
