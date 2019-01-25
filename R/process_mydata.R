#!/usr/bin/env Rscript

# This is my custome script to parse my own data 
# RNA-seq data collection from 4 tomato species

library(dplyr)
source("counts2tpm.R")


# Read S.lycoperiscum read counts dataset
Slyc_data <- read.table("../data/Slyc_mapped_stats.txt", header = TRUE, 
                        sep = "\t", row.names = 1)
Slyc_readCounts <- Slyc_data[,2:ncol(Slyc_data)]
Slyc_TPM <- counts_to_tpm(Slyc_readCounts, Slyc_data$Length)

# prepare vegetative dataset
Slyc_leaf <- rowMeans(Slyc_TPM[ , c("Slyc_leaf_P3m1", "Slyc_leaf_P3m2", "Slyc_leaf_P3m3",
                                    "Slyc_leaf_P4d1", "Slyc_leaf_P4d2", "Slyc_leaf_P4d3",
                                    "Slyc_leaf_P4p1", "Slyc_leaf_P4p2", "Slyc_leaf_P4p3",
                                    "Slyc_leaf_P5d1", "Slyc_leaf_P5d2", "Slyc_leaf_P5d3",
                                    "Slyc_leaf_P6d1", "Slyc_leaf_P6d2", "Slyc_leaf_P6d3",
                                    "Slyc_leaf_P6p1", "Slyc_leaf_P6p2", "Slyc_leaf_P6p3")])
Slyc_vege <- rowMeans(Slyc_TPM[ , c("Slyc_vege1", "Slyc_vege2", "Slyc_vege3")])
Slyc_root <- rowMeans(Slyc_TPM[ , c("Slyc_root1", "Slyc_root2")])
Slyc_stem <- rowMeans(Slyc_TPM[ , c("Slyc_stem1", "Slyc_stem2")])
Slyc_seed <- rowMeans(Slyc_TPM[ , c("Slyc_seed1", "Slyc_seed2")])

Slyc_vegeALL <- data.frame(Slyc_leaf, Slyc_vege, Slyc_root, Slyc_stem, Slyc_seed)
colnames(Slyc_vegeALL) <- c("leaf", "vege_meristem", "root", "stem", "seed")
Slyc_vegeALL <- add_rownames(Slyc_vegeALL, "Genes")

# prepare reproductive dataset
Slyc_flor <- rowMeans(Slyc_TPM[ , c("Slyc_flor1", "Slyc_flor2", "Slyc_flor3")])
Slyc_ovul <- rowMeans(Slyc_TPM[ , c("Slyc_ovul1", "Slyc_ovul2", "Slyc_ovul3", "Slyc_ovul4")])
Slyc_poll <- rowMeans(Slyc_TPM[ , c("Slyc_polm1", "Slyc_polm2", "Slyc_polm3", 
                                    "Slyc_polp1", "Slyc_polp2", "Slyc_polp3")])
Slyc_styl <- rowMeans(Slyc_TPM[ , c("Slyc_stym1", "Slyc_stym2", "Slyc_stym3",
                                    "Slyc_styp1", "Slyc_styp2", "Slyc_styp3")])

Slyc_reprALL <- data.frame(Slyc_flor, Slyc_ovul, Slyc_poll, Slyc_styl)
colnames(Slyc_reprALL) <- c("floral_meristem", "ovule", "pollen", "style")
Slyc_reprALL <- add_rownames(Slyc_reprALL, "Genes")

write.table(Slyc_vegeALL, file = "../output/Slyc_vege_expr.txt", quote = F, row.names = F)
write.table(Slyc_reprALL, file = "../output/Slyc_repr_expr.txt", quote = F, row.names = F)
rm(list = ls(pattern = "Slyc")) 


# Read S.pennellii read counts dataset
Spen_data <- read.table("../data/Spen_mapped_stats.txt", header = TRUE, 
                        sep = "\t", row.names = 1)
Spen_readCounts <- Spen_data[,2:ncol(Spen_data)]
Spen_TPM <- counts_to_tpm(Spen_readCounts, Spen_data$Length)

# prepare vegetative dataset
Spen_leaf <- rowMeans(Spen_TPM[ , c("Spen_leaf_P3m1", "Spen_leaf_P3m2", "Spen_leaf_P3m3",
                                   "Spen_leaf_P4d1", "Spen_leaf_P4d2", "Spen_leaf_P4d3",
                                   "Spen_leaf_P4p1", "Spen_leaf_P4p2", "Spen_leaf_P4p3",
                                   "Spen_leaf_P5d1", "Spen_leaf_P5d2", "Spen_leaf_P5d3",
                                   "Spen_leaf_P6d1", "Spen_leaf_P6d2", "Spen_leaf_P6d3",
                                   "Spen_leaf_P6p1", "Spen_leaf_P6p2", "Spen_leaf_P6p3")])
Spen_root <- rowMeans(Spen_TPM[ , c("Spen_root1", "Spen_root2")])
Spen_stem <- rowMeans(Spen_TPM[ , c("Spen_stem1", "Spen_stem2")])
Spen_seed <- rowMeans(Spen_TPM[ , c("Spen_seed1", "Spen_seed2")])
Spen_vege <- rowMeans(Spen_TPM[ , c("Spen_vege1", "Spen_vege2", "Spen_vege3")])

Spen_vegeALL <- data.frame(Spen_leaf, Spen_root, Spen_stem, Spen_seed, Spen_vege)
colnames(Spen_vegeALL) <- c("leaf", "root", "stem", "seed", "vege_meristem")
Spen_vegeALL <- add_rownames(Spen_vegeALL, "Genes")


# prepare reproductive dataset
Spen_flor <- rowMeans(Spen_TPM[ , c("Spen_flor1", "Spen_flor2")])
Spen_ovul <- rowMeans(Spen_TPM[ , c("Spen_ovul1", "Spen_ovul2", "Spen_ovul3", "Spen_ovul4")])
Spen_poll <- rowMeans(Spen_TPM[ , c("Spen_polm1", "Spen_polm2", "Spen_polm3", 
                                    "Spen_polp1", "Spen_polp2", "Spen_polp3")])
Spen_styl <- rowMeans(Spen_TPM[ , c("Spen_stym1", "Spen_stym2", "Spen_stym3",
                                    "Spen_styp1", "Spen_styp2", "Spen_styp3")])

Spen_reprALL <- data.frame(Spen_flor, Spen_ovul, Spen_poll, Spen_styl)
colnames(Spen_reprALL) <- c("floral_meristem", "ovule", "pollen", "style")
Spen_reprALL <- add_rownames(Spen_reprALL, "Genes")

write.table(Spen_vegeALL, file = "../output/Spen_vege_expr.txt", quote = F, row.names = F)
write.table(Spen_reprALL, file = "../output/Spen_repr_expr.txt", quote = F, row.names = F)
rm(list = ls(pattern = "Spen")) 



# Read S.habrochaites read counts dataset
Shab_data <- read.table("../data/Shab_mapped_stats.txt", header = TRUE, 
                        sep = "\t", row.names = 1)
Shab_readCounts <- Shab_data[,2:ncol(Shab_data)]
Shab_TPM <- counts_to_tpm(Shab_readCounts, Shab_data$Length)

# prepare vegetative dataset
Shab_leaf <- rowMeans(Shab_TPM[ , c("Shab_leaf_P3m1", "Shab_leaf_P3m2", "Shab_leaf_P3m3",
                                    "Shab_leaf_P4d1", "Shab_leaf_P4d2", "Shab_leaf_P4d3",
                                    "Shab_leaf_P4p1", "Shab_leaf_P4p2", "Shab_leaf_P4p3",
                                    "Shab_leaf_P5d1", "Shab_leaf_P5d2", "Shab_leaf_P5d3",
                                    "Shab_leaf_P6d1", "Shab_leaf_P6d2", "Shab_leaf_P6d3",
                                    "Shab_leaf_P6p1", "Shab_leaf_P6p2", "Shab_leaf_P6p3")])
Shab_root <- rowMeans(Shab_TPM[ , c("Shab_root1", "Shab_root2", "Shab_root3",
                               "Shab_root4", "Shab_root5", "Shab_root6")])
Shab_stem <- (Shab_TPM$Shab_stem1)
Shab_seed <- rowMeans(Shab_TPM[ , c("Shab_seed1", "Shab_seed2")])
Shab_vege <- rowMeans(Shab_TPM[ , c("Shab_vege1", "Shab_vege2", "Shab_vege3")])

Shab_vegeALL <- data.frame(Shab_leaf, Shab_root, Shab_stem, Shab_seed, Shab_vege)
colnames(Shab_vegeALL) <- c("leaf", "root", "stem", "seed", "vege_meristem")
Shab_vegeALL <- add_rownames(Shab_vegeALL, "Genes")

# prepare reproductive dataset
Shab_ovul <- (Shab_TPM$Shab_ovul1)
Shab_poll <- rowMeans(Shab_TPM[ , c("Shab_poll1", "Shab_poll2")])
Shab_styl <- rowMeans(Shab_TPM[ , c("Shab_styl1", "Shab_styl2", "Shab_styl3")])

Shab_reprALL <- data.frame(Shab_ovul, Shab_poll, Shab_styl)
colnames(Shab_reprALL) <- c("ovule", "pollen", "style")
Shab_reprALL <- add_rownames(Shab_reprALL, "Genes")

write.table(Shab_vegeALL, file = "../output/Shab_vege_expr.txt", quote = F, row.names = F)
write.table(Shab_reprALL, file = "../output/Shab_repr_expr.txt", quote = F, row.names = F)
rm(list = ls(pattern = "Shab")) 


# Read S.pimpenfolia read counts dataset
Spim_data <- read.table("../data/Spim_mapped_stats.txt", header = TRUE, 
                        sep = "\t", row.names = 1)
Spim_readCounts <- Spim_data[,2:ncol(Spim_data)]
Spim_TPM <- counts_to_tpm(Spim_readCounts, Spim_data$Length)

# prepare vegetative dataset
Spim_leaf <- rowMeans(Spim_TPM[ , c("Spim_leaf_m1", "Spim_leaf_m2", "Spim_leaf_m3", 
                               "Spim_leaf_p1", "Spim_leaf_m2", "Spim_leaf_m3")])
Spim_coty <- rowMeans(Spim_TPM[ , c("Spim_coty1", "Spim_coty2", "Spim_coty3")])
Spim_hypo <- rowMeans(Spim_TPM[ , c("Spim_hypo1", "Spim_hypo2", "Spim_hypo3", "Spim_hypo4")])
Spim_root <- rowMeans(Spim_TPM[ , c("Spim_root1", "Spim_root2", "Spim_root3", "Spim_root4")])
Spim_vege <- rowMeans(Spim_TPM[ , c("Spim_vege1", "Spim_vege2", "Spim_vege3", "Spim_vege4")])
Spim_seed <- Spim_TPM$Spim_seed1

Spim_vegeALL <- data.frame(Spim_leaf, Spim_coty, Spim_hypo, Spim_root, Spim_vege, Spim_seed)
colnames(Shab_vegeALL) <- c("leaf", "coty", "hypo", "root", "vege_meristem", "seed")
Spim_vegeALL <- add_rownames(Spim_vegeALL, "Genes")

# prepare reproductive dataset
Spim_flom <- rowMeans(Spim_TPM[ , c("Spim_flor_m1", "Spim_flor_m2", "Spim_flor_m3", "Spim_flor_m4")])
Spim_flop <- rowMeans(Spim_TPM[ , c("Spim_flor_p1", "Spim_flor_p2", "Spim_flor_p3", "Spim_flor_p4")])
Spim_ovul <- rowMeans(Spim_TPM[ , c("Spim_ovul1", "Spim_ovul2")])

Spim_reprALL <- data.frame(Spim_flom, Spim_flop, Spim_ovul)
colnames(Spim_reprALL) <- c("floral_meristem", "mature_flower", "ovule")
Spim_reprALL <- add_rownames(Spim_reprALL, "Genes")

write.table(Spim_vegeALL, file = "../output/Spim_vege_expr.txt", quote = F, row.names = F)
write.table(Spim_reprALL, file = "../output/Spim_repr_expr.txt", quote = F, row.names = F)
rm(list = ls(pattern = "Spim")) 

