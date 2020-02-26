# pzw
# 20200226

library(phyloseq)
library(pheatmap)
library(ggplot2)
library(plyr)
library(ape)
library(scales)
theme_set(theme_bw())

otumat <- read.csv("OTU_expression.txt", header=TRUE, sep="\t", quote="", row.names=1)
otumat <- as.matrix(otumat)

taxmat <- read.csv("OTU_taxos.txt", header=TRUE, sep="\t", quote="", row.names=1, na.strings="NA")
taxmat <- as.matrix(taxmat)

otu <- otu_table(otumat, taxa_are_rows=TRUE)
tax <- tax_table(taxmat)

physeq <- phyloseq(otu, tax)
physeq_rel <- transform_sample_counts(physeq,function(x)x/sum(x))
TOP30Taxa <- names(sort(taxa_sums(physeq_rel), TRUE)[1:30])
TOP30Taxa <- prune_taxa(TOP30Taxa, physeq_rel)

sampleinfo <- read.csv("sample.txt", sep="\t", header=TRUE, quote="", row.names=1)
sampleinfo <- as.data.frame(sampleinfo)
sampledata <- sample_data(sampleinfo)


random_tree = rtree(ntaxa(physeq), rooted=TRUE, tip.label=taxa_names(physeq))
physeq1 = merge_phyloseq(physeq, sampledata, random_tree)
physeq2 = phyloseq(otu, tax, sampledata, random_tree)

ordu = ordinate(physeq2, "PCoA", "unifrac", weighted=FALSE)
plot_ordination(physeq2, ordu, color="SampleType") + 
  ggtitle("PCoA on weighted-UniFrac distance") + xlab("PC1 [27.5%]") + ylab("PC2 [21.1%]") +
  stat_ellipse(level = 0.95, show.legend = F) + geom_text(aes(label=row.names(sampledata)), size=4)


plot_heatmap(TOP30Taxa, taxa.label="Genus", method = "NMDS", distance = "bray", low = "#990000",
             high = "#006699", na.value="white", trans = log_trans(2))

