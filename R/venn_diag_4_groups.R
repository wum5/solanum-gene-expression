# a function to generate venn diagram of 4 categories of genes
library(VennDiagram)

venn_diag_to_4 <- function(group1, group2, group3, group4){
    x1 <- group1
    x2 <- group2
    x3 <- group3
    x4 <- group4
    x12 <- intersect(group1,group2)
    x13 <- intersect(group1,group3)
    x14 <- intersect(group1,group4)
    x23 <- intersect(group2,group3)
    x24 <- intersect(group2,group4)
    x34 <- intersect(group3,group4)
    x123 <- Reduce(intersect,list(group1,group2,group3))
    x124 <- Reduce(intersect,list(group1,group2,group4))
    x134 <- Reduce(intersect,list(group1,group3,group4))
    x234 <- Reduce(intersect,list(group2,group3,group4))
    x1234 <- Reduce(intersect,list(group1,group2,group3,group4))

    venn.plot <- draw.quad.venn(
      area1 = length(x1),
      area2 = length(x2),
      area3 = length(x3),
      area4 = length(x4),
      n12 = length(x12),
      n13 = length(x13),
      n14 = length(x14),
      n23 = length(x23),
      n24 = length(x24),
      n34 = length(x34),
      n123 = length(x123),
      n124 = length(x124),
      n134 = length(x134),
      n234 = length(x234),
      n1234 = length(x1234),
      category = c("Slyc", "Spen", "Shab", "Spim"),
      fill = c("orange", "red", "green", "blue"),
      lty = "dashed",
      cex = 2,
      cat.cex = 2,
      cat.col = c("orange", "red", "green", "blue")
    )
    
    geneSets <- Reduce(union, list(x123,x234,x124,x134,x1234))
    return(geneSets)
}

