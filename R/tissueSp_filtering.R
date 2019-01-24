
######################################################################## 
### Pull out sexually tissue-specific genes for each species
Slyc_sex_genes <- classify_sex_genes(Slyc_all_tpm, Slyc_repr)
Spen_sex_genes <- classify_sex_genes(Spen_all_tpm, Spen_repr)
Shab_sex_genes <- classify_sex_genes(Shab_all_tpm, Shab_repr)

Slyc_poll_genes <- Slyc_sex_genes[[1]]
Slyc_styl_genes <- Slyc_sex_genes[[2]]
Slyc_ovul_genes <- Slyc_sex_genes[[3]]

Spen_poll_genes <- Spen_sex_genes[[1]]
Spen_styl_genes <- Spen_sex_genes[[2]]
Spen_ovul_genes <- Spen_sex_genes[[3]]

Shab_poll_genes <- Shab_sex_genes[[1]]
Shab_styl_genes <- Shab_sex_genes[[2]]
Shab_ovul_genes <- Shab_sex_genes[[3]]

## pollen-specific (show specific expression in at least two out three species)
poll12 <- intersect(Slyc_poll_genes,Spen_poll_genes)
poll13 <- intersect(Slyc_poll_genes,Shab_poll_genes)
poll23 <- intersect(Spen_poll_genes,Shab_poll_genes)
poll123 <- Reduce(intersect, list(Slyc_poll_genes,Spen_poll_genes,Shab_poll_genes))
poll_specific <- Reduce(union, list(poll12,poll13,poll23,poll123))

## style-specific (show specific expression in at least two out three species)
styl12 <- intersect(Slyc_styl_genes,Spen_styl_genes)
styl13 <- intersect(Slyc_styl_genes,Shab_styl_genes)
styl23 <- intersect(Spen_styl_genes,Shab_styl_genes)
styl123 <- Reduce(intersect, list(Slyc_styl_genes,Spen_styl_genes,Shab_styl_genes))
styl_specific <- Reduce(union, list(styl12,styl13,styl23,styl123))

## pollen-specific (show specific expression in at least two out three species)
ovul12 <- intersect(Slyc_ovul_genes,Spen_ovul_genes)
ovul13 <- intersect(Slyc_ovul_genes,Shab_ovul_genes)
ovul23 <- intersect(Spen_ovul_genes,Shab_ovul_genes)
ovul123 <- Reduce(intersect, list(Slyc_ovul_genes,Spen_ovul_genes,Shab_ovul_genes))
ovul_specific <- Reduce(union, list(ovul12,ovul13,ovul23,ovul123))

length(poll_specific)
length(styl_specific)
length(ovul_specific)

write.table(poll_specific, file = "pollen_specific_genes.txt", quote=FALSE, row.names=FALSE)
write.table(styl_specific, file = "style_specific_genes.txt", quote=FALSE, row.names=FALSE)
write.table(ovul_specific, file = "ovule_specific_genes.txt", quote=FALSE, row.names=FALSE)
######################################################################## 


