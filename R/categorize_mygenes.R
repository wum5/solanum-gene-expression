#!/usr/bin/env Rscript

# This script is to classify genes into three categories:
# reproductive genes, vegetative genes or generalized genes
source("genes_classifier.R")
source("genes_intersecter.R")


# Pull out vege/repr-specific genes for each species
Slyc_vege_tpm <- read.table("../output/Slyc_vege_expr.txt", header = TRUE, row.names = 1)
Slyc_repr_tpm <- read.table("../output/Slyc_repr_expr.txt", header = TRUE, row.names = 1)
Slyc_classified_genes <- classify_genes(Slyc_vege_tpm, Slyc_repr_tpm)
Slyc_repr <- Slyc_classified_genes[[1]]
Slyc_vege <- Slyc_classified_genes[[2]]
Slyc_genr <- Slyc_classified_genes[[3]]
# number of reproductive genes in Slyc
length(Slyc_repr)
# number of vegetative genes in Slyc
length(Slyc_vege)
# number of generalized genes in Slyc
length(Slyc_genr)


Spen_vege_tpm <- read.table("../output/Spen_vege_expr.txt", header = TRUE, row.names = 1)
Spen_repr_tpm <- read.table("../output/Spen_repr_expr.txt", header = TRUE, row.names = 1)
Spen_classified_genes <- classify_genes(Spen_vege_tpm, Spen_repr_tpm)
Spen_repr <- Spen_classified_genes[[1]]
Spen_vege <- Spen_classified_genes[[2]]
Spen_genr <- Spen_classified_genes[[3]]
# number of reproductive genes in Spen
length(Spen_repr)
# number of vegetative genes in Spen
length(Spen_vege)
# number of generalized genes in Spen
length(Spen_genr)


Spim_vege_tpm <- read.table("../output/Spim_vege_expr.txt", header = TRUE, row.names = 1)
Spim_repr_tpm <- read.table("../output/Spim_repr_expr.txt", header = TRUE, row.names = 1)
Spim_classified_genes <- classify_genes(Spim_vege_tpm, Spim_repr_tpm)
Spim_repr <- Spim_classified_genes[[1]]
Spim_vege <- Spim_classified_genes[[2]]
Spim_genr <- Spim_classified_genes[[3]]
# number of reproductive genes in Spim
length(Spim_repr)
# number of vegetative genes in Spim
length(Spim_vege)
# number of generalized genes in Spim
length(Spim_genr)


Shab_vege_tpm <- read.table("../output/Shab_vege_expr.txt", header = TRUE, row.names = 1)
Shab_repr_tpm <- read.table("../output/Shab_repr_expr.txt", header = TRUE, row.names = 1)
Shab_classified_genes <- classify_genes(Shab_vege_tpm, Shab_repr_tpm)
Shab_repr <- Shab_classified_genes[[1]]
Shab_vege <- Shab_classified_genes[[2]]
Shab_genr <- Shab_classified_genes[[3]]
# number of reproductive genes in Shab
length(Shab_repr)
# number of vegetative genes in Shab
length(Shab_vege)
# number of generalized genes in Shab
length(Shab_genr)


# generate the Venn diagram showing reproductive-specific genes shared across 4 species
# also return the genes shared by at least 3 species
pdf("../output/repr_genes_venn.pdf", width=5, height=4)
reprGenes <- venn_diag_3by4(Slyc_repr, Spen_repr, Shab_repr, Spim_repr)
dev.off()

length(reprGenes)
write.table(reprGenes, file = "../output/shared_RPgenes.txt", sep="\t", 
            quote = F, row.names = F, col.names = F)


# generate the Venn diagram showing vegetative-specific genes shared across 4 species
# also return the genes shared by at least 3 species
pdf("../output/vege_genes_venn.pdf", width=5, height=4)
vegeGenes <- venn_diag_3by4(Slyc_vege, Spen_vege, Shab_vege, Spim_vege)
dev.off()

length(vegeGenes)
write.table(vegeGenes, file = "../output/shared_VEgenes.txt", sep="\t", 
            quote = F, row.names = F, col.names = F)


# generate the Venn diagram showing generalized genes shared across 4 species
# also return the genes shared by at least 3 species
pdf("../output/genr_genes_venn.pdf", width=5, height=4)
genrGenes <- venn_diag_3by4(Slyc_genr, Spen_genr, Shab_genr, Spim_genr)
dev.off()

length(genrGenes)
write.table(genrGenes, file = "../output/shared_GRgenes.txt", sep="\t", 
            quote = F, row.names = F, col.names = F)



# Pull out sexually tissue-specific genes for each species
Slyc_sex_genes <- classify_sex_genes(Slyc_repr_tpm, Slyc_repr)
Spen_sex_genes <- classify_sex_genes(Spen_repr_tpm, Spen_repr)
Shab_sex_genes <- classify_sex_genes(Shab_repr_tpm, Shab_repr)

Slyc_poll_genes <- Slyc_sex_genes[[1]]
Slyc_styl_genes <- Slyc_sex_genes[[2]]
Slyc_ovul_genes <- Slyc_sex_genes[[3]]

Spen_poll_genes <- Spen_sex_genes[[1]]
Spen_styl_genes <- Spen_sex_genes[[2]]
Spen_ovul_genes <- Spen_sex_genes[[3]]

Shab_poll_genes <- Shab_sex_genes[[1]]
Shab_styl_genes <- Shab_sex_genes[[2]]
Shab_ovul_genes <- Shab_sex_genes[[3]]

# pollen-specific (show specific expression in at least two out three species)
poll_specf <- intersct_2by3(Slyc_poll_genes, Spen_poll_genes, Shab_poll_genes)
length(poll_specf)

# style-specific (show specific expression in at least two out three species)
styl_specf <- intersct_2by3(Slyc_styl_genes, Spen_styl_genes, Shab_styl_genes)
length(styl_specf)

# pollen-specific (show specific expression in at least two out three species)
ovul_specf <- intersct_2by3(Slyc_ovul_genes, Spen_ovul_genes, Shab_ovul_genes)
length(ovul_specf)

# write the output to tables
write.table(poll_specf, file = "../output/pollen_specific_genes.txt", 
            quote = F, row.names = F, col.names = F)
write.table(styl_specf, file = "../output/style_specific_genes.txt", 
            quote = F, row.names = F, col.names = F)
write.table(ovul_specf, file = "../output/ovule_specific_genes.txt", 
            quote = F, row.names = F, col.names = F)


