library("phyloseq")
library("vegan")
library("DESeq2")
library("ggplot2")
library("dendextend")
library("tidyr")
library("viridis")
library("reshape")

setwd("D:/practice/R/plotAMP")
otumat <- read.csv("D:/practice/R/plotAMP/ASV/otumat.csv", header=TRUE, quote="", row.names=1)
taxmat <- read.csv("D:/practice/R/plotAMP/ASV/taxmat.csv", header=TRUE, quote="", row.names=1)

otu <- otu_table(as.matrix(otumat), taxa_are_rows=FALSE)
tax <- tax_table(as.matrix(taxmat))
physeq = phyloseq(otu, tax)

# diversity bar plot
plot_bar(physeq, fill="Phylum")

# alpha
plot_richness(physeq, measures=c("Chao1", "Shannon"))










