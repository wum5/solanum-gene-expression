# A function to use TPM<0.5 and TPM>2 as cutoffs filtering out for genes to be
# 1) reproductive-specific; 2) vegetative-specific; and 3) generalized.

classify_genes <- function(vege_tissues, repr_tissues) {
  keep <- rowSums(vege_tissues>0.5) == 0
  vege_keep <- rownames(vege_tissues[keep, ])
  keep <- rowSums(repr_tissues>2) >= 1
  repr_keep <- rownames(repr_tissues[keep, ])
  repr_genes <- intersect(vege_keep, repr_keep)
  length(repr_genes)
  
  keep <- rowSums(vege_tissues>2) >= 1
  vege_keep <- rownames(vege_tissues[keep, ])
  keep <- rowSums(repr_tissues>0.5) == 0
  repr_keep <- rownames(repr_tissues[keep, ])
  vege_genes <- intersect(vege_keep, repr_keep)
  length(vege_genes)
  
  keep <- rowSums(vege_tissues>2) == length(head(vege_tissues,1))
  vege_keep <- rownames(vege_tissues[keep, ])
  keep <- rowSums(repr_tissues>2) == length(head(repr_tissues,1))
  repr_keep <- rownames(repr_tissues[keep, ])
  genr_genes <- intersect(vege_keep, repr_keep)
  length(genr_genes)
  
  return(list(repr_genes, vege_genes, genr_genes))
}


# A function to filter out for male/female-specific genes, but have a option
# to allow they are expressed in some number (cutoff here) of tissues 
classify_sex_genes <- function(tpm_table, repr_genes){
  if ("polm" %in% colnames(tpm_table))
  { pollen <- data.frame(tpm_table$polm, tpm_table$polp) }
  else { pollen <- data.frame(tpm_table$poll, tpm_table$poll) }
  
  if ("stym" %in% colnames(tpm_table))
  { style <- data.frame(tpm_table$stym,tpm_table$styp) }
  else { style <- data.frame(tpm_table$styl, tpm_table$styl) }
  
  ovule <- data.frame(tpm_table$ovul, tpm_table$ovul)
  
  rownames(pollen) <- rownames(tpm_table)  
  rownames(style) <- rownames(tpm_table)  
  rownames(ovule) <- rownames(tpm_table)  
  
  # pollen-specific
  keep <- rowSums(ovule>0.5) == 0
  filter1 <- rownames(pollen[keep, ])
  keep <- rowSums(style>0.5) == 0
  filter2 <- rownames(pollen[keep, ])  
  keep <- rowSums(pollen>2) >= 1
  filter3 <- rownames(pollen[keep, ])
  poll_genes <- Reduce(intersect,list(filter1, filter2, filter3, repr_genes))
  length(poll_genes)
  
  # style-specific
  keep <- rowSums(ovule>0.5) == 0
  filter1 <- rownames(pollen[keep, ])
  keep <- rowSums(pollen>0.5) == 0
  filter2 <- rownames(pollen[keep, ])  
  keep <- rowSums(style>2) >= 1
  filter3 <- rownames(style[keep, ])
  styl_genes <- Reduce(intersect,list(filter1, filter2, filter3, repr_genes))
  length(poll_genes)
  
  # ovule-specific 
  keep <- rowSums(pollen>0.5) == 0
  filter1 <- rownames(pollen[keep, ])
  keep <- rowSums(style>0.5) == 0
  filter2 <- rownames(pollen[keep, ])  
  keep <- rowSums(ovule>2) >= 1
  filter3 <- rownames(ovule[keep, ])
  ovul_genes <- Reduce(intersect,list(filter1, filter2, filter3, repr_genes))
  length(poll_genes)
  
  return(list(poll_genes, styl_genes, ovul_genes))
}

