library("phyloseq")
library("vegan")
library("DESeq2")
library("ggplot2")
library("dendextend")
library("tidyr")
library("viridis")
library("reshape")

setwd("D:/practice/R/plotAMP")
count_tab <- read.csv("D:/practice/R/plotAMP/ASV/otumat.csv", header=T, row.names=1, check.names=F)
tax_tab <- as.matrix(read.table("D:/practice/R/plotAMP/ASV/taxmat.csv", header=T, row.names=1, check.names=F, na.strings="", sep=","))
sample_info_tab <- read.table("D:/practice/R/plotAMP/ASV/sample_info.txt", header=T, row.names=1, check.names=F)

# Alpha
rarecurve(t(count_tab), step=100, col=sample_info_tab$color, lwd=2, ylab="ASVs")
abline(v=(min(rowSums(t(filt_count_tab)))))





















