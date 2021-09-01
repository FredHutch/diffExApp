# create a sample matrix from input condition names and selected column names -----
data <- read.table("~/Documents/work/dataCore/shiny/diffEx/data/counts.csv", sep = ",", header = TRUE)
condition1_name <- "wt"
condition2_name <- "trt"
condition1_selected <- c("SRR1039508", "SRR1039509", "SRR1039512")
condition2_selected <- c("SRR1039521", "SRR1039520", "SRR1039517")
gene_col <- "gene"

###################
## DESEQ HELPERS ##
###################
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
  idx <- rowSums( counts(dds, normalized=TRUE) >= 5 ) >= 3
  dds <- dds[idx]
  
  # run deseq analysis
  deseq <- DESeq(dds)
  
  return(deseq)
}

###################
## MAIN FUNCTION ##
###################

deseq <- function(data,
                  condition1_name,
                  condition2_name,
                  condition1_selected,
                  condition2_selected,
                  gene_col) {
  
  # Create the sample matrix
  sample_matrix <- createSampleMatrix(condition1_name,
                                      condition2_name,
                                      condition1_selected,
                                      condition2_selected)
  
  # prep data so that it matches sample_matrix
  prepped_data <- prepDataset(data,
                              sample_matrix,
                              gene_col)
  
  # run deseq analysis
  # output is deseq data obj
  de_out <- deseqDE(prepped_data,
                    sample_matrix)
  
  # get results from de_out
  de_res <- data.frame(results(de_out))
  
  return(de_res)
  
}
