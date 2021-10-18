## GENERIC DE HELPERS #############################################################################
# create a sample matrix from selected column names and supplied naming info
createSampleMatrix <- function(condition1_name,
                               condition2_name,
                               condition1_selected,
                               condition2_selected) {
  samples <- c(condition1_selected, condition2_selected)
  condition <- c(rep(condition1_name, length(condition1_selected)),
                 rep(condition2_name, length(condition2_selected)))
  
  # create sample matrix
  # FIXME: This will error out if a replicate is selected twice, should I create a cleaner error catch?
  # Try to break app with this
  sample_matrix <- data.frame(row.names = samples, condition = condition)
  
  # deseq expects factors, cond1 always first
  sample_matrix$condition <- factor(sample_matrix$condition, levels = c(condition1_name, condition2_name))
  
  return(sample_matrix)
}

# subset dataset and makes sure it matches sample matrix
prepDataset <- function(data,
                        sample_matrix,
                        gene_col) {
  # make gene_col row names of data
  data <- column_to_rownames(data, gene_col)
  
  # check that all samples in sample_matrix are in colnames of data
  if (!all(row.names(sample_matrix) %in% colnames(data))) {
    # if they aren't all present, error out
    stop("ERROR: sample_matrix contains samples that are not present in counts column names")
  }
  
  # check that all colnames in data are in sample matrix
  if (!all(colnames(data) %in% row.names(sample_matrix))) {
    # if not remove extra columns from data
    data <- data[,colnames(data) %in% row.names(sample_matrix)]
  }
  
  # check that colnames of data are in the same order as sample matrix
  if (!all(rownames(sample_matrix) == colnames(data))) {
    # if not reorder columns in data
    data <- data[rownames(sample_matrix)]
  }
  
  return(data)
}

## DESEQ DE #######################################################################################
# run deseq de analysis
# output deseq dataset object
deseqDE <- function(data,
                    sample_matrix) {
  # using sample matrix and data create deseq dataset object
  dds <- DESeqDataSetFromMatrix(countData = data,
                                colData = sample_matrix,
                                design = ~condition)
  
  # filter out low counts
  dds <- estimateSizeFactors(dds)
  idx <- rowSums(counts(dds, normalized=TRUE) >= 5) >= 3
  dds <- dds[idx]
  
  # run deseq analysis
  deseq <- DESeq(dds)
  
  return(deseq)
}

## EDGER DE #######################################################################################
edgerDE <- function(data,
                    sample_matrix) {
  # create dgelist
  dge <- DGEList(counts = data, group = sample_matrix$condition)
  
  # filter out genes with low normalized counts
  # keeping metrics same as deseq (at least three columns > 5)
  idx <- rowSums(cpm(dge) >= 5) >= 3
  dge$counts <- dge$counts[idx, ]
  
  # recalculate library size after filtering low counts
  dge$samples$lib.size <- colSums(dge$counts)
  
  # TMM normalization is applied to this dataset to account for compositional difference between
  dge <- calcNormFactors(dge)
  
  # Before we fit negative binomial GLMs, we need to define our design matrix based on the experimental design
  condition <- sample_matrix$condition
  design <- model.matrix(~condition)
  rownames(design) <- rownames(sample_matrix)
  
  # estimate the NB dispersion
  dge <- estimateDisp(dge, design, robust=TRUE)
  
  # Determine differentially expressed genes
  # Fit genewise glms
  fit <- glmFit(dge, design)
  
  # Conduct likelihood ratio test
  lrt <- glmLRT(fit)
  
  # include norm counts in output
  lrt$counts <- dge$counts
  lrt$normalizedCounts <- cpm(dge$counts)
  
  return(lrt)
}

## MAIN DE FUNCTION ###############################################################################
diffEx <- function(data,
                   condition1_name,
                   condition2_name,
                   condition1_selected,
                   condition2_selected,
                   gene_col,
                   de_package) {
  
  # create sample matrix
  sample_matrix <- createSampleMatrix(condition1_name,
                                      condition2_name,
                                      condition1_selected,
                                      condition2_selected) 
  # prep dataset
  prepped_data <- prepDataset(data,
                              sample_matrix,
                              gene_col)
  
  # run de analysis
  if (de_package == "DESeq2") {
    de_out <- deseqDE(prepped_data, sample_matrix)
  } else {
    de_out <- edgerDE(prepped_data, sample_matrix)
  }
  
  return(de_out)
}

## FILTER / FORMAT DE_OUT #########################################################################
# get counts based on user settings for de_package
getCounts <- function(de_out,
                      de_package,
                      normalized = TRUE) {
  if  (de_package == "DESeq2") {
    counts <- counts(de_out, normalized = normalized)
    } else if (de_package == "edgeR") {
      if (normalized) {
      counts <- data.frame(de_out$normalizedCounts)
      } else {
      counts <- data.frame(de_out$counts)
      }}
  
  return(counts)
}


# filter differential expression output based on a set pvalue and logfc threshold
getResults <- function(de_out,
                       de_package) {
  # get results from de output
  if (de_package == "DESeq2") {
    de_res <- data.frame(results(de_out))
  } else if (de_package == "edgeR") {
    de_res <- data.frame(topTags(de_out, n = Inf))
  }
  
  return(de_res)
}

# format results based on user inputs
formatResults <- function(de_res,
                          logfc_threshold,
                          pvalue_threshold,
                          logfc_col,
                          pvalue_col,
                          de_column,
                          de_filter,
                          subset) {
  
  # create index of de genes
  de_idx <- abs(de_res[[logfc_col]]) >= logfc_threshold & de_res[[pvalue_col]] <= pvalue_threshold
  
  # DESEQ for some reasons turns pval/padj to NA, this messes up the de_idx (NA instead of T/F)
  # https://bioconductor.org/packages/release/bioc/vignettes/DESeq2/inst/doc/DESeq2.html#pvaluesNA
  # change NA to FALSE in de_idx so these cases are dropped
  de_idx[is.na(de_idx)] = FALSE
  
  if (de_filter == "both") {
    # save col in results
    de_res$isDE <- de_idx
  } else if (de_filter == "up") {
    # find up values
    up_idx <- de_res[[logfc_col]] >= 0
    # calc where both vectors == TRUE
    up_de_idx <- rowSums(data.frame(de_idx, up_idx)) == 2
    # save col in results
    de_res$isDE <- up_de_idx
  } else if (de_filter == "down") {
    # find up values
    dn_idx <- de_res[[logfc_col]] <= 0
    # calc where both vectors == TRUE
    dn_de_idx <- rowSums(data.frame(de_idx, dn_idx)) == 2
    # save col in results
    de_res$isDE <- dn_de_idx
    }

  if (subset) {
    de_res <- de_res[de_res$isDE, ]
  }
  
  # if de_column = FALSE, remove isDE column 
  if (!de_column) {
    de_res <- subset(de_res, select= -c(isDE))
  }
  
  return(de_res)
  
}
