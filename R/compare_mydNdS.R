#!/usr/bin/env Rscript

# This is my custome script to plot the dN/dS distribution of genes 
# in different categories based on their expression profiles

library(ggplot2)

# read PAML output (i.e. dN/dS data)
dndsData <- read.table("../data/PAML_out_4sp_rt1.txt", sep = "\t", header=TRUE)
dndsData <- dndsData[dndsData[2]<5,]
head(dndsData)
nrow(dndsData)

# read categorized genes
RPgenes <- read.table("../output/shared_RPgenes.txt")
VEgenes <- read.table("../output/shared_VEgenes.txt")
GRgenes <- read.table("../output/shared_GRgenes.txt")

Ovlgenes <- read.table("../output/ovule_specific_genes.txt")
Polgenes <- read.table("../output/pollen_specific_genes.txt")
Stygenes <- read.table("../output/style_specific_genes.txt")

# convert the input format for following processing
RP_list <- cbind(sub("*\\.[0-9]", "", RPgenes$V1), "RP")
VE_list <- cbind(sub("*\\.[0-9]", "", VEgenes$V1), "VE")
GR_list <- cbind(sub("*\\.[0-9]", "", GRgenes$V1), "GR")

Ovl_list <- cbind(sub("*\\.[0-9]", "", Ovlgenes$V1), "Ovl")
Pol_list <- cbind(sub("*\\.[0-9]", "", Polgenes$V1), "Pol")
Sty_list <- cbind(sub("*\\.[0-9]", "", Stygenes$V1), "Sty")

all_genes <- data.frame(rbind(RP_list, VE_list, GR_list))
colnames(all_genes) <- c("Genes", "Group")
head(all_genes)

repr_genes <- data.frame(rbind(Ovl_list, Pol_list, Sty_list))
colnames(repr_genes) <- c("Genes", "Group")
head(repr_genes)

# combine two dataframes by gene ids
all_dNdS <- merge(all_genes, dndsData, by="Genes")
# number of reproductive-specific genes
nrow(all_genes[all_genes$Group=="RP",] )
# number of vegetative-specific genes
nrow(all_genes[all_genes$Group=="VE",] )
# number of generalized genes
nrow(all_genes[all_genes$Group=="GR",] )

repr_dNdS <- merge(repr_genes, dndsData, by="Genes")
# number of ovule-specific genes
nrow(repr_genes[repr_dNdS$Group=="Ovl",] )
# number of pollen-specific genes
nrow(repr_genes[repr_dNdS$Group=="Pol",] )
# number of style-specific genes
nrow(repr_genes[repr_dNdS$Group=="Sty",] )

# plot the dN/dS distribution of all genes in 3 categories 
pdf("../output/dNdS_allGenes.pdf", width=6, height=4)
ggplot(all_dNdS, aes(dNdS, fill = Group)) + 
  geom_density(alpha = 0.2) +
  scale_fill_manual(values=c("skyblue", "hotpink", "green")) +
  scale_x_continuous(limits = c(0, 3)) +
  ggtitle("dN/dS distribution of genes from 3 categories")
dev.off()

tapply(all_dNdS$dNdS, all_dNdS$Group, summary)

# plot the dN/dS distribution of reproductive-specific genes in 3 tissues 
pdf("../output/dNdS_reprGenes.pdf", width=6, height=4)
ggplot(repr_dNdS, aes(dNdS, fill = Group)) + 
  geom_density(alpha = 0.2) +
  scale_fill_manual(values=c("skyblue", "hotpink", "green")) +
  scale_x_continuous(limits = c(0, 5)) +
  ggtitle("dN/dS distribution of tissue-specific genes")
dev.off()

tapply(repr_dNdS$dNdS, repr_dNdS$Group, summary)


