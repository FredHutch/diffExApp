# plot de output in MA plot
# basemean on x axis
# log2foldchange on y axis
resultsToMa <- function(results,
                        de_package,
                        pvalue_threshold,
                        logfc_threshold,
                        fdr) {
  
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
  de_idx <- abs(results[[logfc_col]]) >= logfc_threshold & results[[pvalue_col]] <= pvalue_threshold
  # NA to FALSE
  de_idx[is.na(de_idx)] = FALSE
  # Reset levels so TRUE is first
  de_idx <- factor(de_idx, levels = c(TRUE, FALSE))
  
  # plot
  ma <- ggplot(results, aes(.data[[mean_counts_col]], .data[[logfc_col]], colour = de_idx)) +
    geom_point(size=2, alpha = .7) + 
    scale_y_continuous(limits=c(-3, 3), oob = squish)
  
  # if package is DESeq log10 the base mean axis
  if (de_package == "DESeq2") {
    ma <- ma +
      scale_x_log10()
  }
  
  ma +
    geom_hline(yintercept = 0, colour="grey", size=1) + 
    labs(x="mean of normalized counts", y="log fold change") + 
    theme_classic()
}