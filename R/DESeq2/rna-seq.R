library(DESeq2)
library(ggplot2)
library(pheatmap)
library(dplyr)


# 读入
data <- read.table("final_featureCounts.txt", header=TRUE, skip=1, row.names=1)
colnames(data) <- gsub(".bam", "", colnames(data), fixed=TRUE)
colnames(data) <- gsub("bam.", "", colnames(data), fixed=TRUE)
countdata <- data[ , 6:ncol(data)]

# 计算TPM
KB <- data$Length / 1000
RPK <- countdata / KB
TPM <- t(t(RPK) / colSums(RPK) * 1000000)
TPM <- merge(data, as.data.frame(TPM), by="row.names", sort=FALSE)
TPM <- TPM[, 1:ncol(TPM)]
write.table(TPM, "final_featureCounts.TPM.txt", sep="\t", quote=FALSE, row.names=FALSE)

# 读入metadata
metadata <- read.table("metadata.txt", header=TRUE)
rownames(metadata) <- metadata$sampleid
metadata[match(colnames(countdata), metadata$sampleid), ]

# 计算差异
dds <- DESeqDataSetFromMatrix(countData=countdata, colData=metadata, design=~Group)
dds <- dds[ rowSums(counts(dds)) > 1, ]
dds <- DESeq(dds)
res <- results(dds, pAdjustMethod="fdr", alpha=0.05)
res <- res[order(res$padj),]
summary(res)
mcols(res, use.names=TRUE)

# 保存结果
resdata <- merge(as.data.frame(res), as.data.frame(counts(dds, normalized=TRUE)),
                 by="row.names", sort=FALSE)
write.csv(resdata, file="LoGlu_HiGlu_mm39_diff.csv", row.names=FALSE)


# 也可以输出差异基因
diff_gene <- subset(res, padj < 0.05 & (log2FoldChange > 1 | log2FoldChange < -1))
write.table(x=as.data.frame(diff_gene), file="results_gene_annotated_significant.txt", sep="\t", quote=F, col.names=NA)


# PCA
# Convert all samples to rlog
dds_rlog <- rlog(dds, blind=FALSE)

# Plot PCA by column variable
plotPCA(dds_rlog, intgroup="Group", ntop=500) +
  theme_bw() + # remove default ggplot2 theme
  geom_point(size=5) + # Increase point size
  scale_y_continuous(limits=c(-5, 5)) + # change limits to fix figure dimensions
  ggtitle(label="Principal Component Analysis (PCA)",
          subtitle="Top 500 most variable genes")



# heatmap
# Convert all samples to rlog
dds_rlog <- rlog(dds, blind=FALSE)

# Gather 30 significant genes and make matrix
mat <- assay(dds_rlog[row.names(diff_gene)])[1:40, ]

# Choose which column variables you want to annotate the columns by.
annotation_col <- data.frame(
  Group=factor(colData(dds_rlog)$Group),
  Replicate=factor(colData(dds_rlog)$Replicate),
  row.names=colData(dds_rlog)$sampleid
)

# Specify colors you want to annotate the columns by.
ann_colors <- list(
  Group=c(LoGlu="lightblue", HiGlu="darkorange"),
  Replicate=c(Rep1="darkred", Rep2="forestgreen")
)

# Make Heatmap with pheatmap function.
## See more in documentation for customization
pheatmap(mat=mat, 
         color=colorRampPalette(brewer.pal(9, "YlOrBr"))(255), 
         scale="row", # Scale genes to Z-score (how many standard deviations)
         annotation_col=annotation_col, # Add multiple annotations to the samples
         annotation_colors=ann_colors,# Change the default colors of the annotations
         fontsize=6.5, # Make fonts smaller
         cellwidth=55, # Make the cells wider
         show_colnames=F)

# 火山图
vol_data <- data.frame(gene=row.names(res), pval=-log10(res$padj), lfc=res$log2FoldChange)

# Remove any rows that have NA as an entry
vol_data <- na.omit(vol_data)

# Color the points which are up or down
## If fold-change > 0 and pvalue > 1.3 (Increased significant)
## If fold-change < 0 and pvalue > 1.3 (Decreased significant)
vol_data <- mutate(vol_data, color = case_when(vol_data$lfc > 0 & vol_data$pval > 1.3 ~ "Increased",
                                               vol_data$lfc < 0 & vol_data$pval > 1.3 ~ "Decreased",
                                               vol_data$pval < 1.3 ~ "nonsignificant"))

# Make a basic ggplot2 object with x-y values
vol <- ggplot(vol_data, aes(x=lfc, y=pval, color=color))

# Add ggplot2 layers
vol +   
  ggtitle(label="Volcano Plot", subtitle="Colored by fold-change direction") +
  geom_point(size=2.5, alpha=0.8, na.rm=T) +
  scale_color_manual(name="Directionality",
                     values=c(Increased="#008B00", Decreased="#CD4F39", nonsignificant="darkgray")) +
  theme_bw(base_size=14) + # change overall theme
  theme(legend.position="right") + # change the legend
  xlab(expression(log[2]("LoGlu" / "HiGlu"))) + # Change X-Axis label
  ylab(expression(-log[10]("adjusted p-value"))) + # Change Y-Axis label
  geom_hline(yintercept=1.3, colour="darkgrey") + # Add p-adj value cutoff line
  scale_y_continuous(trans="log1p") # Scale yaxis due to large p-values


# MA
plotMA(res, ylim=c(-5, 5))

## Dispersions
plotDispEsts(dds)

## 单基因
# Convert all samples to rlog
dds_rlog <- rlog(dds, blind=FALSE)

# Get gene with highest expression
top_gene <- rownames(res)[which.min(res$log2FoldChange)]

# Plot single gene
plotCounts(dds=dds, 
           gene=top_gene, 
           intgroup="Group", 
           normalized=T, 
           transform=T)






