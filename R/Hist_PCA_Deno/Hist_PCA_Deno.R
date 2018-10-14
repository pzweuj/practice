# 20181014

brain12 <- read.table("Allen_BrainAtlas_12regions_Microarray.txt", header=TRUE)

### Histograms
# 先设置一下12种颜色
col <- c("red",
         "blue",
         "orange",
         "yellow",
         "magenta",
         "khaki",
         "green",
         "gray",
         "goldenrod",
         "cyan",
         "cornsilk",
         "coral")

# 应该可以利用随机函数还有colors()函数来随机生成12种颜色，更方便

# 设定输出目标
pdf("Histograms.pdf", width=8.27, height=11.69)

# 合并下面的图
par(mfrow=c(4, 3))

# 循环出图
for(i in c(1:12)){
  xname = colnames(brain12[i])
  hist(brain12[, xname], main=paste("Histogram of ", xname),
       col=col[i],
       xlab="Probe Intensity", ylab="Frequency")
  rm(xname)
}

# 输出结束
dev.off()

### PCA PLOTS
# install.packages("ggfortify")
library(ggfortify)

# 前处理
brain12 <- t(brain12)
brain12 <- as.data.frame(brain12)
brain12$region <- rownames(brain12)
brain12[is.na(brain12)] <- 0

# 计算pca
brain12.data <- brain12[c(1:2817)]
pca <- prcomp(brain12.data)

summary <- summary(pca)

# 输出pdf
pdf("PCA.pdf")
autoplot(pca, data=brain12, colour="region")
dev.off()

# 输出summary
write.table(summary$importance,
            file="PCA_summary.txt",
            sep='\t', col.names=NA, quote=FALSE)


### HIERARCHICAL CLUSTERING
# Scale the data
brain12.data.scale <- scale(brain12.data)

# Use Pearson’s correlations
brain12.data.cor <- cor(t(brain12.data.scale), method="pearson")

# calculate distances
brain12.data.dist <- dist(brain12.data.cor)

# Use the “single” method as the agglomeration (linkage) method
hc <- hclust(brain12.data.dist, method="single")

pdf("Dendrogram.pdf")
plot(hc, main="Cluster Dendrogram")
dev.off()
