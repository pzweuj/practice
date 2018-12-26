library(ggfortify)
library(factoextra)

###
GSE57820 <- read.table("GSE57820_non_normalized.txt", header=T, sep="\t", row.names=1)
samples <- seq(from=1, to=23, by=2)
GSE57820 <- GSE57820[, samples]
namelist <- list()
for(i in colnames(GSE57820)){
  j <- unlist(strsplit(i, "_"))
  k <- paste(j[2], j[3], j[4], sep="_")
  namelist <- append(namelist, k)
}
colnames(GSE57820) <- unlist(namelist)
GSE57820.df <- as.data.frame(t(GSE57820))
GSE57820.data <- scale(GSE57820.df)
GSE57820.df$samples <- rownames(GSE57820.df)

## PCA
pca <- prcomp(GSE57820.data)
autoplot(pca, data=GSE57820.df, colour="samples")

## HC
dist.cor <- get_dist(GSE57820.data, method="pearson")
# 这个出来k=1，不知道咋办
fviz_nbclust(GSE57820.data, kmeans, method="gap_stat")
hc.complete.cor <- hclust(d=dist.cor, method="complete")
fviz_dend(hc.complete.cor, labels=T)



