library(ChIPseeker)
library(TxDb.Hsapiens.UCSC.hg19.knownGene)
txdb <- TxDb.Hsapiens.UCSC.hg19.knownGene
library(clusterProfiler)
library(org.Hs.eg.db)
library(ReactomePA)

peak <- readPeakFile("Peak/H3k27ac_peaks.narrowPeak", header=FALSE)
peak

# peaks在基因组中的位置
covplot(peak, weightCol="V5")
covplot(peak, weightCol="V5", chrs=c("chr17", "chr18"), xlim=c(4.5e7, 5e7))

# 与TSS区域结合的peaks的概况
promoter <- getPromoters(TxDb=txdb, upstream=3000, downstream=3000)
tagMatrix <- getTagMatrix(peak, windows=promoter)

# ChIP与TSS区域结合的热图
tagHeatmap(tagMatrix, xlim=c(-3000, 3000), color="red")

# ChIP peak的平均谱与TSS区域结合
plotAvgProf(tagMatrix, xlim=c(-3000, 3000),
            xlab="Genomic Region (5'->3')",
            ylab = "Read Count Frequency")

plotAvgProf(tagMatrix, xlim=c(-3000, 3000), conf = 0.95, resample = 1000)

## anno
peakAnno <- annotatePeak(peak, tssRegion=c(-3000, 3000),
                         TxDb=txdb, annoDb="org.Hs.eg.db")

write.table(as.data.frame(peakAnno), "Anno/peaks.anno.xls", quote=F, row.names=F, sep="\t")

# 饼图
plotAnnoPie(peakAnno)

# 条形图
plotAnnoBar(peakAnno)

# 韦恩pie
vennpie(peakAnno)

# upset
upsetplot(peakAnno)

# combine
upsetplot(peakAnno, vennpie=TRUE)

# 可视化TF结合基因座相对于TSS的分布
plotDistToTSS(peakAnno,
              title="Distribution of transcription factor-binding loci\nrelative to TSS")


# 富集分析
pathway1 <- enrichPathway(as.data.frame(peakAnno)$geneId)
head(pathway1, 2)
gene <- seq2gene(peak, tssRegion = c(-1000, 1000), flankDistance = 3000, TxDb=txdb)
pathway2 <- enrichPathway(gene)
head(pathway2, 2)
dotplot(pathway2)

#
q("no")