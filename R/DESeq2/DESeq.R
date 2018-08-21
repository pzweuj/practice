library(DESeq2)
library(ggplot2)

# 数据预处理
sampleNames <- c("MCF10A_1", "MCF10A_2", "MCF10A_3", "MCF7_1", "MCF7_2", "MCF7_3")
data <- read.table("all_feature.txt", header=T, quote="\t", skip=1)
names(data)[7:12] <- sampleNames
countData <- as.matrix(data[7:12])
rownames(countData) <- data$Geneid
database <- data.frame(name=sampleNames, condition=c("MCF10A", "MCF10A", "MCF10A", "MCF7", "MCF7", "MCF7"))
rownames(database) <- sampleNames

# 设置分组信息并构建dds对象
dds <- DESeqDataSetFromMatrix(countData, colData=database, design= ~ condition)
dds <- dds[ rowSums(counts(dds)) > 1, ]

# 使用DESeq函数估计离散度，然后差异分析获得res对象
dds <- DESeq(dds)
res <- results(dds)
write.csv(res, "res_output.csv")
resdata <- merge(as.data.frame(res), as.data.frame(counts(dds, normalized=TRUE)),by="row.names",sort=FALSE)
write.csv(resdata, file="all_output.csv", row.names=F)

# PCA
rld <- rlog(dds)
pcaData <- plotPCA(rld, intgroup=c("condition", "name"), returnData=T)
percentVar <- round(100*attr(pcaData, "percentVar"))
pca <- ggplot(pcaData, aes(PC1, PC2, color=condition, shape=name)) + 
	geom_point(size=3) + 
	ggtitle("DESeq2 PCA") + 
	xlab(paste0("PC1: ", percentVar[1], "% variance")) + 
	ylab(paste0("PC2: ", percentVar[2], "% variance"))

# MA
library(geneplotter)
plotMA(res, main="DESeq2", ylim=c(-2, 2))

# valcano
resdata$change <- as.factor(
	ifelse(
		resdata$padj<0.01 & abs(resdata$log2FoldChange)>1,
		ifelse(resdata$log2FoldChange>1, "Up", "Down"),
		"NoDiff"
	)
)

valcano <- ggplot(data=resdata, aes(x=log2FoldChange, y=-log10(padj), color=change)) + 
	geom_point(alpha=0.8, size=1) + 
	theme_bw(base_size=15) + 
	theme(
		panel.grid.minor=element_blank(),
		panel.grid.major=element_blank()
	) + 
	ggtitle("DESeq2 Valcano") + 
	scale_color_manual(name="", values=c("red", "green", "black"), limits=c("Up", "Down", "NoDiff")) + 
	geom_vline(xintercept=c(-1, 1), lty=2, col="gray", lwd=0.5) + 
	geom_hline(yintercept=-log10(0.01), lty=2, col="gray", lwd=0.5)

# heatmap
library(pheatmap)

sum(res$padj<0.1, na.rm=T)
select <- order(rowMeans(counts(dds, normalized=T)), decreasing=T)[1:1000]
nt <- normTransform(dds)
log2.norm.counts <- assay(nt)[select,]
df <- as.data.frame(colData(dds)[, c("name", "condition")])
pheatmap(log2.norm.counts, cluster_rows=T, show_rownames=F, cluster_cols=T, annotation_col=df, fontsize=6)
