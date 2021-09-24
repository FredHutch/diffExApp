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
resultsToMa <- function(de_res,
                        logfc_col,
                        mean_counts_col,
                        x_label,
                        y_label,
                        legend_title,
                        de_vec) {
  
  # Reorder de_vec factors so that TRUE is first
  de_vec <- factor(de_vec, levels=c("TRUE", "FALSE"))
  
  # base plot
  ma <- ggplot(de_res, aes(.data[[mean_counts_col]], .data[[logfc_col]])) +
    geom_point(size=2, alpha = .6, aes(colour = de_vec)) + 
    scale_y_continuous(limits=c(-3, 3), oob = squish)
  
  # if package is DESeq log10 the base mean axis
  # check for basemean col
  basemean_present <- "baseMean" %in% names(de_res)
  if (basemean_present) {
    ma <- ma +
      scale_x_log10()
  }
  
  # add hline
  ma <- ma +
    geom_hline(yintercept = 0, colour="grey", size=1, linetype = "dashed")
  
  # add custom axes and legend labels
  ma <- ma +
    labs(x = x_label, y = y_label, color = legend_title) +
    theme_classic(base_size = 12)
  
  # display plot
  ma
}
