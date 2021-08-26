# create a sample matrix from input condition names and selected column names -----
# condition1_name <- "wt"
# condition2_name <- "trt"
# select_condition1 <- c("SRR1039508", "SRR1039509", "SRR1039512")
# select_condition2 <- c("SRR1039521", "SRR1039520", "SRR1039517")

createSampleMatrix <- function(condition1_name,
                               condition2_name,
                               select_condition1,
                               select_condition2) {
  samples <- c(select_condition1, select_condition2)
  condition <- c(rep(condition1_name, length(select_condition1)),
                 rep(condition2_name, length(select_condition2)))
  sample_mat <- data.frame(row.names = samples, condition = condition)
  return(sample_mat)
}