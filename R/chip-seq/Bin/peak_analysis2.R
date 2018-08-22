library(ChIPseeker)
library(TxDb.Mmusculus.UCSC.mm10.knownGene)
txdb <- TxDb.Mmusculus.UCSC.mm10.knownGene
library(clusterProfiler)
library(org.Mm.eg.db)
library(ReactomePA)

cbx7 <- readPeakFile("Peak/cbx7_peaks.narrowPeak", header=FALSE)
Ring1B <- readPeakFile("Peak/Ring1B_peaks.narrowPeak", header=FALSE)
SUZ12 <- readPeakFile("Peak/SUZ12_peaks.narrowPeak", header=FALSE)

peak_list <- list(cbx7=cbx7, Ring1B=Ring1B, SUZ12=SUZ12)

# peaks在基因组中的位置
covplot(cbx7, weightCol="V5")
covplot(Ring1B, weightCol="V5")
covplot(SUZ12, weightCol="V5")
covplot(cbx7, weightCol="V5", chrs=c("chr17", "chr18"), xlim=c(4.5e7, 5e7))

# 与TSS区域结合的peaks的概况
promoter <- getPromoters(TxDb=txdb, upstream=3000, downstream=3000)
tagMatrix <- lapply(peak_list, getTagMatrix, windows=promoter)

# ChIP与TSS区域结合的热图
tagHeatmap(tagMatrix, xlim=c(-3000, 3000), color=NULL)

# ChIP peak的平均谱与TSS区域结合
plotAvgProf(tagMatrix, xlim=c(-3000, 3000),
            xlab="Genomic Region (5'->3')",
            ylab = "Read Count Frequency", facet="row")

# plotAvgProf(tagMatrix, xlim=c(-3000, 3000), conf = 0.95, resample = 1000)

## anno
cbx7_Anno <- annotatePeak(cbx7, tssRegion=c(-3000, 3000),
                         TxDb=txdb, annoDb="org.Mm.eg.db", verbose=FALSE)

Ring1B_Anno <- annotatePeak(Ring1B, tssRegion=c(-3000, 3000),
                         TxDb=txdb, annoDb="org.Mm.eg.db", verbose=FALSE)

SUZ12_Anno <- annotatePeak(SUZ12, tssRegion=c(-3000, 3000),
                         TxDb=txdb, annoDb="org.Mm.eg.db", verbose=FALSE)

write.table(as.data.frame(cbx7_Anno), "Anno/cbx7.anno.xls", quote=F, row.names=F, sep="\t")
write.table(as.data.frame(Ring1B_Anno), "Anno/Ring1B.anno.xls", quote=F, row.names=F, sep="\t")
write.table(as.data.frame(SUZ12_Anno), "Anno/SUZ12.anno.xls", quote=F, row.names=F, sep="\t")

anno_list <- lapply(peak_list, annotatePeak, TxDb=txdb, tssRegion=c(-3000, 3000), verbose=FALSE, annoDb="org.Mm.eg.db")

# single sample
plotAnnoPie(cbx7_Anno)
vennpie(cbx7_Anno)
upsetplot(cbx7_Anno)
upsetplot(cbx7_Anno, vennpie=TRUE)

# 条形图
plotAnnoBar(anno_list)

# 可视化TF结合基因座相对于TSS的分布
plotDistToTSS(anno_list,
              title="Distribution of transcription factor-binding loci\nrelative to TSS")


# 富集分析
genes = lapply(anno_list, function(i) as.data.frame(i)$geneId)
names(genes) = sub("_", "\n", names(genes))
compKEGG <- compareCluster(geneCluster=genes, fun="enrichKEGG", pvalueCutoff=0.05, pAdjustMethod="BH")
dotplot(compKEGG, showCategory=15, title="KEGG Pathway Enrichment Analysis")

#
genes = lapply(anno_list, function(i) as.data.frame(i)$geneId)
vennplot(genes)

#
q("no")
