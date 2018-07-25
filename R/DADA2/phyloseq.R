library(phyloseq)
library(ggplot2)

theme_set(theme_bw())

otumat <- read.csv("~/workspace/AMP/OTU/otumax.nochim.csv", header=TRUE, quote="", row.names=1)
taxmat <- read.csv("~/workspace/AMP/OTU/taxmat.csv", header=TRUE, quote="", row.names=1)

otumat <- as.matrix(otumat)
taxmat <- as.matrix(taxmat)

otu <- otu_table(otumat, taxa_are_rows=FALSE)
tax <- tax_table(taxmat)

physeq = phyloseq(otu, tax)

# bar plot
plot_bar(physeq, fill="Domain")
plot_bar(physeq, fill="Phylum")
plot_bar(physeq, fill="Class")
plot_bar(physeq, fill="Order")
plot_bar(physeq, fill="Family")
plot_bar(physeq, fill="Genus")
plot_bar(physeq, fill="Species")

# alpha
plot_richness(physeq, measures="Shannon")
plot_richness(physeq, measures="Simpson")
plot_richness(physeq, measures=c("Shannon", "Simpson"))

# beta
names <- rownames(otu)
physeq.pop <- transform_sample_counts(physeq, function(otu) otu/sum(otu))
ord.nmds.bray <- ordinate(physeq.pop, method="NMDS", distance="bray")
plot_ordination(physeq.pop, ord.nmds.bray, title="Bray NMDS")
