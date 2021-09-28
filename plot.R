## COUNTS TABLE VISUALIZATIONS ####################################################################
# plot PCA
plotPca <- function(counts,
                    sample_matrix) {
  
  # compute PCA on transposed counts
  pca <- prcomp(t(counts))
  
  # pull out info on PCs as dataframe
  pca_df <- data.frame(pca$x)
  
  # add col for annotations
  pca_df$condition <- as.character(sample_matrix$condition)
  
  # plot
  ggplot(pca_df, aes(x = PC1, y = PC2, color = condition)) +
    geom_point() +
    ggrepel::geom_label_repel(aes(label = rownames(pca_df)), show.legend = FALSE) +
    theme_classic()
}

# plot heatmap
plotHeatmap <- function(counts,
                        sample_matrix,
                        de_vec,
                        silent = TRUE) {
  # subset de genes only
  de_counts <- counts[de_vec, ]
  
  # heatmap
  pheatmap(de_counts, 
           scale= "row", 
           border_color = NA, 
           annotation_col = sample_matrix,
           silent = silent)
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
