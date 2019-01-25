# A function to use cutoffs (default: TPM1_cut=0.5 and TPM2_cut=2) filtering out for genes 
# to be 1) reproductive-specific; 2) vegetative-specific; or 3) generalized.

# the minimum cutoff for defining expressed genes
TPM1_cut <- 0.5 
# the maximum cutoff for defining non-expressed genes
TPM2_cut <- 2


classify_genes <- function(vege_tissues, repr_tissues) {
    keep <- rowSums(vege_tissues>TPM1_cut) == 0
    vege_keep <- rownames(vege_tissues[keep, ])
    keep <- rowSums(repr_tissues>TPM2_cut) >= 1
    repr_keep <- rownames(repr_tissues[keep, ])
    repr_genes <- intersect(vege_keep, repr_keep)
    length(repr_genes)
  
    keep <- rowSums(vege_tissues>TPM2_cut) >= 1
    vege_keep <- rownames(vege_tissues[keep, ])
    keep <- rowSums(repr_tissues>TPM1_cut) == 0
    repr_keep <- rownames(repr_tissues[keep, ])
    vege_genes <- intersect(vege_keep, repr_keep)
    length(vege_genes)
  
    keep <- rowSums(vege_tissues>TPM2_cut) == length(head(vege_tissues,1))
    vege_keep <- rownames(vege_tissues[keep, ])
    keep <- rowSums(repr_tissues>TPM2_cut) == length(head(repr_tissues,1))
    repr_keep <- rownames(repr_tissues[keep, ])
    genr_genes <- intersect(vege_keep, repr_keep)
    length(genr_genes)
  
    return(list(repr_genes, vege_genes, genr_genes))
}


# A function to filter out for male/female-specific genes, but have a option
# to allow they are expressed in some number (cutoff here) of tissues 
classify_sex_genes <- function(tpm_table, repr_genes){
    pollen <- tpm_table[,"pollen", drop = F]
    style <- tpm_table[,"style", drop = F]
    ovule <- tpm_table[,"ovule", drop = F]
  
    rownames(pollen) <- rownames(tpm_table)  
    rownames(style) <- rownames(tpm_table)  
    rownames(ovule) <- rownames(tpm_table)  
  
    # pollen-specific
    keep <- rowSums(ovule>TPM1_cut) == 0
    filter1 <- rownames(pollen[keep, , drop=F])
    keep <- rowSums(style>TPM1_cut) == 0
    filter2 <- rownames(pollen[keep, , drop=F])  
    keep <- rowSums(pollen>TPM2_cut) >= 1
    filter3 <- rownames(pollen[keep, , drop=F])
    poll_genes <- Reduce(intersect,list(filter1, filter2, filter3, repr_genes))
    length(poll_genes)
  
    # style-specific
    keep <- rowSums(ovule>TPM1_cut) == 0
    filter1 <- rownames(pollen[keep, , drop=F])
    keep <- rowSums(pollen>TPM1_cut) == 0
    filter2 <- rownames(pollen[keep, , drop=F])  
    keep <- rowSums(style>TPM2_cut) >= 1
    filter3 <- rownames(style[keep, , drop=F])
    styl_genes <- Reduce(intersect,list(filter1, filter2, filter3, repr_genes))
    length(poll_genes)
  
    # ovule-specific 
    keep <- rowSums(pollen>TPM1_cut) == 0
    filter1 <- rownames(pollen[keep, , drop=F])
    keep <- rowSums(style>TPM1_cut) == 0
    filter2 <- rownames(pollen[keep, , drop=F])  
    keep <- rowSums(ovule>TPM2_cut) >= 1
    filter3 <- rownames(ovule[keep, , drop=F])
    ovul_genes <- Reduce(intersect,list(filter1, filter2, filter3, repr_genes))
    length(poll_genes)
  
    return(list(poll_genes, styl_genes, ovul_genes))
}


