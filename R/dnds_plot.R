setwd("/Users/mengwu/Documents/Research/Solanum_RP")
library(ggplot2)

dndsData <- read.table("dndsValues.txt", sep = "\t", header=TRUE)
dndsData <- dndsData[dndsData[2]<5,]
head(dndsData)

#and combine into your new data frame vegLengths
ggplot(dndsData, aes(dNdS, fill = tissue)) + 
  geom_density(alpha = 0.2) +
  scale_fill_manual(values=c("skyblue", "hotpink", "green")) +
  scale_x_continuous(limits = c(0, 2))

tapply(dndsData$dNdS, dndsData$tissue, summary)

nrow(RPgenes)
RPgenes <- dndsData[dndsData$tissue=='RP',]
nrow(VEgenes)
VEgenes <- dndsData[dndsData$tissue=='VE',]
nrow(GRgenes)
GRgenes <- dndsData[dndsData$tissue=='GR',]


dndsData <- read.table("dndsValues-2.txt", sep = "\t", header=TRUE)
dndsData <- dndsData[dndsData[2]<5,]
head(dndsData)

ggplot(dndsData, aes(x=dNdS, fill=tissue)) + xlim(0,3) +
  geom_histogram(alpha=0.5, position="identity", colour="black", bins = 30) +
  theme_classic()



