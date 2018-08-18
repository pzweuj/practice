library("phyloseq")
library("vegan")
library("DESeq2")
library("ggplot2")
library("dendextend")
#library("tidyr")
#library("viridis")
#library("reshape")

countdata <- read.table("ASV/ASVs_counts.txt", header=T, row.names=1, check.names=F)
taxdata <- as.matrix(read.table("ASV/ASVs_taxonomy.txt", header=T, row.names=1, check.names=F, na.strings="", sep="\t"))
sample_info <- read.table("Rawdata/sample_info.txt", header=T, row.names=1, check.names=F)
# taxdata <- replace(taxdata, taxdata == " ", NA)

## 数据初始化
# 排序
# ord_col <- as.numeric(sapply(names(countdata), FUN=function(x){strsplit(x, "F3D")[[1]][2]}))
# countdata <- countdata[order(ord_col)]

# ord_row <- as.numeric(sapply(row.names(sample_info), FUN=function(x){strsplit(x, "F3D")[[1]][2]}))
# sample_info <- sample_info[order(ord_row),]

# 设置颜色
sample_info$color[sample_info$time == "Early"] <- "green"
sample_info$color[sample_info$time == "Late"] <- "red"

# 创建新的phyloseq对象
count_phy <- otu_table(countdata, taxa_are_rows=T)
tax_phy <- tax_table(taxdata)
sample_info_phy <- sample_data(sample_info)
ASV_physeq <- phyloseq(count_phy, tax_phy, sample_info_phy)


## 分类分析
plot_bar(ASV_physeq, fill="Phylum") + 
  theme_bw() +
  theme(axis.text.x=element_text(angle=90, vjust=0.4, hjust=1), legend.title=element_blank()) +
  labs(x="Sample", y="Copies recovered", title="All samples")

## Alpha分析
# 稀释曲线
rarecurve(t(countdata), step=100, col=sample_info$color, lwd=2, ylab="# of ASVs", xlab="# of Sequences")
abline(v=(min(rowSums(t(countdata)))))

# 画出Chao1和Shannon
# 横坐标为样本
plot_richness(ASV_physeq, color="time", measures=c("Chao1", "Shannon")) + 
  scale_color_manual(values=unique(sample_info$color[order(sample_info$time)])) +
  theme(legend.title = element_blank())
# 横坐标为状态
plot_richness(ASV_physeq, x="dpw", color="time", measures=c("Chao1", "Shannon")) + 
  scale_color_manual(values=unique(sample_info$color[order(sample_info$time)])) +
  theme(legend.title = element_blank())


## Beta分析
# 使用DESeq2进行差异分析
deseq_counts <- DESeqDataSetFromMatrix(countdata, colData = sample_info, design = ~time)
deseq_counts_vst <- varianceStabilizingTransformation(deseq_counts)

# 生成表格
vst_trans_count <- assay(deseq_counts_vst)

# 计算出差异矩阵
dist_count <- dist(t(vst_trans_count))

# 聚类
clust_count <- hclust(dist_count, method="ward.D2")
dend_count <- as.dendrogram(clust_count, hang=0.1)
dend_cols <- sample_info$color[order.dendrogram(dend_count)]
labels_colors(dend_count) <- dend_cols

# 画出树
plot(dend_count, ylab="VST count. dist.")

# 使用差异矩阵
vst_count_phy <- otu_table(vst_trans_count, taxa_are_rows=T)
vst_physeq <- phyloseq(vst_count_phy, sample_info_phy)

# 使用phyloseq画PCoA图
vst_pcoa <- ordinate(vst_physeq, method="MDS", distance="euclidean")
eigen_vals <- vst_pcoa$values$Eigenvalues

# PCoA图
plot_ordination(vst_physeq, vst_pcoa, color="time") + 
  labs(col="dpw") + geom_point(size=1) + 
  geom_text(aes(label=rownames(sample_info), hjust=0.3, vjust=-0.4)) + 
  coord_fixed(sqrt(eigen_vals[2]/eigen_vals[1])) + ggtitle("PCoA") + 
  scale_color_manual(values=unique(sample_info$color[order(sample_info$time)])) + 
  theme(legend.position="none")


##
quit("no")








