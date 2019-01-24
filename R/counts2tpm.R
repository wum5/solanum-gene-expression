# A function to calculates TPM from raw reads counts by dividing the raw read counts 
# by the gene size and then by the summed norlimized total counts (per million)

counts_to_tpm <- function(counts, geneLength) {
  
  # Ensure valid arguments.
  stopifnot(length(geneLength) == nrow(counts))
  
  # Process one column at a time.
  tpm <- do.call(cbind, lapply(1:ncol(counts), function(i) {
    rate = log(counts[,i]) - log(geneLength[i])
    denom = log(sum(exp(rate)))
    exp(rate - denom + log(1e6))
  }))
  
  # Copy the row and column names from the original matrix.
  colnames(tpm) <- colnames(counts)
  rownames(tpm) <- rownames(counts)  
  tpm <- data.frame(tpm)
  return(tpm)
}


