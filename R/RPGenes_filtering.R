setwd("/Users/mengwu/Documents/Research/Solanum_RP")
source("genes_classify_functions.R")
source("venn_diag_4_groups.R")


########################################################################
### Pull out vege/repr-specific genes for each species
Slyc_vege_tpm <- read.table("Slyc_vege_expr.txt", header = TRUE, row.names = 1)
Slyc_repr_tpm <- read.table("Slyc_repr_expr.txt", header = TRUE, row.names = 1)
Slyc_all_tpm  <- cbind(Slyc_vege_tpm, Slyc_repr_tpm)
Slyc_classfied_genes <- classify_genes(Slyc_vege_tpm, Slyc_repr_tpm)
Slyc_repr <- Slyc_classfied_genes[[1]]
Slyc_vege <- Slyc_classfied_genes[[2]]
Slyc_genr <- Slyc_classfied_genes[[3]]
  
Spen_vege_tpm <- read.table("Spen_vege_expr.txt", header = TRUE, row.names = 1)
Spen_repr_tpm <- read.table("Spen_repr_expr.txt", header = TRUE, row.names = 1)
Spen_all_tpm  <- cbind(Spen_vege_tpm, Spen_repr_tpm)
Spen_classfied_genes <- classify_genes(Spen_vege_tpm, Spen_repr_tpm)
Spen_repr <- Spen_classfied_genes[[1]]
Spen_vege <- Spen_classfied_genes[[2]]
Spen_genr <- Spen_classfied_genes[[3]]

Spim_vege_tpm <- read.table("Spim_vege_expr.txt", header = TRUE, row.names = 1)
Spim_repr_tpm <- read.table("Spim_repr_expr.txt", header = TRUE, row.names = 1)
Spim_all_tpm  <- cbind(Spim_vege_tpm, Spim_repr_tpm)
Spim_classfied_genes <- classify_genes(Spim_vege_tpm, Spim_repr_tpm)
Spim_repr <- Spim_classfied_genes[[1]]
Spim_vege <- Spim_classfied_genes[[2]]
Spim_genr <- Spim_classfied_genes[[3]]

Shab_vege_tpm <- read.table("Shab_vege_expr.txt", header = TRUE, row.names = 1)
Shab_repr_tpm <- read.table("Shab_repr_expr.txt", header = TRUE, row.names = 1)
Shab_all_tpm  <- cbind(Shab_vege_tpm, Shab_repr_tpm)
Shab_classfied_genes <- classify_genes(Shab_vege_tpm, Shab_repr_tpm)
Shab_repr <- Shab_classfied_genes[[1]]
Shab_vege <- Shab_classfied_genes[[2]]
Shab_genr <- Shab_classfied_genes[[3]]

### Venn diagram of reproductive-specific genes
reprGenes <- venn_diag_to_4(Slyc_repr, Spen_repr, Shab_repr, Spim_repr)

### generate the output
write.table(reprGenes, file = "shared_RPgenes.txt", sep="\t", 
            quote=FALSE, row.names=FALSE)


### Venn diagram of vege-specific genes
vegeGenes <- venn_diag_to_4(Slyc_vege, Spen_vege, Shab_vege, Spim_vege)

### generate the output
write.table(vegeGenes, file = "shared_VGgenes.txt", sep="\t", 
            quote=FALSE, row.names=FALSE)


# Venn diagram of generalized-expressed genes
genrGenes <- venn_diag_to_4(Slyc_genr, Spen_genr, Shab_genr, Spim_genr)

# generate the output
write.table(genrGenes, file = "shared_GRgenes.txt", sep="\t", 
            quote=FALSE, row.names=FALSE)




