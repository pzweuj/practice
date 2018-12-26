library(ggfortify)
library(factoextra)

# 数据处理
GSE57820 <- read.table("GSE57820_series_matrix.txt",
                       skip=72, row.names=1, header=T, nrow=47323)
GSE57820 <- as.data.frame(t(GSE57820))
info <- read.table("info.txt", header=T)
GSE57820$samples <- info$title
GSE57820.data <- GSE57820[c(1:47323)]

# PCA
pca <- prcomp(GSE57820.data)
autoplot(pca, data=GSE57820, colour="samples")

# HC
rownames(GSE57820.data) <- GSE57820$samples
dist.cor <- get_dist(scale(GSE57820.data), method="pearson")
# k=1，不知道咋办
fviz_nbclust(GSE57820.data, kmeans, method="gap_stat")
hc.complete.cor <- hclust(d=dist.cor, method="complete")
fviz_dend(hc.complete.cor, labels=T)