library("phyloseq")
library("vegan")
library("DESeq2")
library("ggplot2")
library("dendextend")
library("tidyr")
library("viridis")
library("reshape")

countdata <- read.table("DADA2/ASVs_counts.txt", header=T, row.names=1, check.names=F)
taxdata <- as.matrix(read.table("DADA2/ASVs_taxonomy.txt", header=T, row.names=1, check.names=F, na.strings="", sep="\t"))

sample_info <- read.table("Rawdata/sample_info.txt", header=T, row.names=1, check.names=F)

## 把样本分一下类，之前说过了是有空白对照的
# 把对照组和实验组的总count数统计出来，方便之后做差异分析
blank_ASV_counts <- rowSums(countdata[,1:4])
sample_ASV_counts <- rowSums(countdata[,5:20])

# 由于实验组是16个样本，对照组是4个样本，所以把实验组的结果除以4来标准化
norm_sample_ASV_counts <- sample_ASV_counts/4

# 找出大概的差异基因
blank_ASVs <- names(blank_ASV_counts[blank_ASV_counts * 10 > norm_sample_ASV_counts])
length(blank_ASVs)

# 去掉空白之后剩下的
colSums(countdata[!rownames(countdata) %in% blank_ASVs, ]) / colSums(countdata) * 100

# 利用刚刚的结果过滤一次
filt_count_tab <- countdata[!rownames(countdata) %in% blank_ASVs, -c(1:4)]
# 过滤之后剩下的样本的信息
filt_sample_info_tab<-sample_info[-c(1:4), ]

# 设置作图颜色
filt_sample_info_tab$color[filt_sample_info_tab$char == "water"] <- "blue"
filt_sample_info_tab$color[filt_sample_info_tab$char == "biofilm"] <- "darkgreen"
filt_sample_info_tab$color[filt_sample_info_tab$char == "altered"] <- "chocolate4"
filt_sample_info_tab$color[filt_sample_info_tab$char == "glassy"] <- "black"
filt_sample_info_tab$color[filt_sample_info_tab$char == "carbonate"] <- "darkkhaki"

filt_sample_info_tab

### BETA分析
# 使用DESeq2进行差异分析
deseq_counts <- DESeqDataSetFromMatrix(filt_count_tab, colData = filt_sample_info_tab, design = ~type)
deseq_counts_vst <- varianceStabilizingTransformation(deseq_counts)

# 生成表格
vst_trans_count_tab <- assay(deseq_counts_vst)

# 计算出差异矩阵
euc_dist <- dist(t(vst_trans_count_tab))

# 聚类
euc_clust <- hclust(euc_dist, method="ward.D2")
euc_dend <- as.dendrogram(euc_clust, hang=0.1)
dend_cols <- filt_sample_info_tab$color[order.dendrogram(euc_dend)]
labels_colors(euc_dend) <- dend_cols

# 画出树
plot(euc_dend, ylab="VST Euc. dist.")

# 使用差异矩阵
vst_count_phy <- otu_table(vst_trans_count_tab, taxa_are_rows=T)
sample_info_tab_phy <- sample_data(filt_sample_info_tab)
vst_physeq <- phyloseq(vst_count_phy, sample_info_tab_phy)

# 使用phyloseq画PCoA图
vst_pcoa <- ordinate(vst_physeq, method="MDS", distance="euclidean")
eigen_vals <- vst_pcoa$values$Eigenvalues

# PCoA图
plot_ordination(vst_physeq, vst_pcoa, color="char") + 
  labs(col="type") + geom_point(size=1) + 
  geom_text(aes(label=rownames(filt_sample_info_tab), hjust=0.3, vjust=-0.4)) + 
  coord_fixed(sqrt(eigen_vals[2]/eigen_vals[1])) + ggtitle("PCoA") + 
  scale_color_manual(values=unique(filt_sample_info_tab$color[order(filt_sample_info_tab$char)])) + 
  theme(legend.position="none")

### ALPHA分析
# 稀释曲线
rarecurve(t(filt_count_tab), step=100, col=filt_sample_info_tab$color, lwd=2, ylab="# of ASVs", xlab="# of Sequences")
abline(v=(min(rowSums(t(filt_count_tab)))))

# 富集还有丰度分析
# 创建新的phyloseq对象
count_tab_phy <- otu_table(filt_count_tab, taxa_are_rows=T)
tax_tab_phy <- tax_table(taxdata)

ASV_physeq <- phyloseq(count_tab_phy, tax_tab_phy, sample_info_tab_phy)

# 画出Chao1和Shannon
# 横坐标为样本
plot_richness(ASV_physeq, color="char", measures=c("Chao1", "Shannon")) + 
  scale_color_manual(values=unique(filt_sample_info_tab$color[order(filt_sample_info_tab$char)])) +
  theme(legend.title = element_blank())
# 横坐标为状态
plot_richness(ASV_physeq, x="type", color="char", measures=c("Chao1", "Shannon")) + 
  scale_color_manual(values=unique(filt_sample_info_tab$color[order(filt_sample_info_tab$char)])) +
  theme(legend.title = element_blank())

###其他分类分析
# 门分类
phyla_counts_tab <- otu_table(tax_glom(ASV_physeq, taxrank="Phylum")) 
phyla_tax_vec <- as.data.frame(tax_table(tax_glom(ASV_physeq, taxrank="Phylum")))$Phylum
rownames(phyla_counts_tab) <- as.vector(phyla_tax_vec)

# we also have to account for sequences that weren't assigned any taxonomy even at the phylum level 
# these came into R as 'NAs' in the taxonomy table, but their counts are still in the count table
# so we can get that value for each sample by substracting the column sums of this new table (that has everything that had a phylum assigned to it) from the column sums of the starting count table (that has all representative sequences)
unannotated_tax_counts <- colSums(filt_count_tab) - colSums(phyla_counts_tab)
# and we'll add this row to our phylum count table:
phyla_and_unidentified_counts_tab <- rbind(phyla_counts_tab, "Unannotated"=unannotated_tax_counts)

# now we'll remove the Proteobacteria, so we can next add them back in broken down by class
temp_major_taxa_counts_tab <- phyla_and_unidentified_counts_tab[!row.names(phyla_and_unidentified_counts_tab) %in% "Proteobacteria", ]

class_counts_tab <- otu_table(tax_glom(ASV_physeq, taxrank="Class")) # making count table broken down by class (contains classes beyond the Proteobacteria too at this point)
class_tax_tab <- data.frame(tax_table(tax_glom(ASV_physeq, taxrank="Class"))) # getting a table of these class names
proteo_classes_vec <- as.vector(class_tax_tab$Class[class_tax_tab$Phylum == "Proteobacteria"]) # making a vector of just the Proteobacteria classes

rownames(class_counts_tab) <- as.vector(class_tax_tab$Class) # changing the row names like above so that they correspond to the taxonomy, rather than an ASV identifier
proteo_class_counts_tab <- class_counts_tab[row.names(class_counts_tab) %in% proteo_classes_vec, ]

# there are also likely some some sequences that were resolved to the level of Proteobacteria, but not any further, and therefore would be missing from our class table
# we can find the sum of them by subtracting the proteo clas count table from just the Proteobacteria row from the original phylum-level count table
proteo_no_class_annotated_counts <- phyla_and_unidentified_counts_tab[row.names(phyla_and_unidentified_counts_tab) %in% "Proteobacteria", ] - colSums(proteo_class_counts_tab)

# now combining the tables:
major_taxa_counts_tab <- rbind(temp_major_taxa_counts_tab, proteo_class_counts_tab, "Unresolved_Proteobacteria"=proteo_no_class_annotated_counts)
# and to check we didn't miss any other sequences, we can compare the column sums to see if they are the same:
identical(colSums(major_taxa_counts_tab), colSums(filt_count_tab)) # if "TRUE", we know nothing fell through the cracks

# now we'll generate a proportions table for summarizing:
major_taxa_proportions_tab <- apply(major_taxa_counts_tab, 2, function(x) x/sum(x)*100)

# if we check the dimensions of this table at this point
dim(major_taxa_proportions_tab)
# we see there are currently 25 rows, which might be a little busy for a summary figure
# many of these taxa make up a very small percentage, so we're going to filter some out
# this is a completely arbitrary decision solely to ease visualization and intepretation, entirely up to your data and you
# here, we'll only keep rows (taxa) that make up greater than 5% in any individual sample
temp_filt_major_taxa_proportions_tab <- data.frame(major_taxa_proportions_tab[apply(major_taxa_proportions_tab, 1, max) > 5, ])
# checking how many we have that were above this threshold
dim(temp_filt_major_taxa_proportions_tab) # now we have 13, much more manageable for an overview figure

# though each of the filtered taxa made up less than 5% alone, together they can be more
# so we're going to add a row called "Other" that keeps track of how much we filtered out (which will also keep our totals at 100%)
filtered_proportions <- colSums(major_taxa_proportions_tab) - colSums(temp_filt_major_taxa_proportions_tab)
filt_major_taxa_proportions_tab <- rbind(temp_filt_major_taxa_proportions_tab, "Other"=filtered_proportions)

