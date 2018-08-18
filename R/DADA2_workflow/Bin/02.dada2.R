library(dada2) # loads DADA2

list.files('Filtdata')

# setting a few variables we're going to use
fnFs <- sort(list.files('Filtdata', pattern="_sub_R1_trim_.fq.gz", full.names=TRUE))
fnRs <- sort(list.files('Filtdata', pattern="_sub_R2_trim_.fq.gz", full.names=TRUE))

sample.names <- sapply(strsplit(basename(fnFs), "_"), `[`, 1)

filtFs <- file.path("Cleandata", paste0(sample.names, "_sub_R1_filtered.fq.gz"))
filtRs <- file.path("Cleandata", paste0(sample.names, "_sub_R2_filtered.fq.gz"))

plotQualityProfile(fnFs)
plotQualityProfile(fnRs)
plotQualityProfile(fnFs[17:20])

out <- filterAndTrim(
  fnFs, filtFs, fnRs, filtRs, truncLen=c(250, 200),
  maxN=0, maxEE=c(2, 2), truncQ=2, rm.phix=TRUE, minLen=175,
  compress=TRUE, multithread=TRUE # 在windows下，multithread设置成FALSE
)
head(out)

class(out) # matrix
dim(out) # 20 2

plotQualityProfile(filtFs)
plotQualityProfile(filtRs)
plotQualityProfile(filtFs[17:20])

# 分别计算正向和反向
errF <- learnErrors(filtFs, multithread=TRUE)
errR <- learnErrors(filtRs, multithread=TRUE)
# 画出错误率统计图
plotErrors(errF, nominalQ=TRUE)
plotErrors(errR, nominalQ=TRUE)

derepFs <- derepFastq(filtFs, verbose=TRUE)
derepRs <- derepFastq(filtRs, verbose=TRUE)
# 用sample.names来改名
names(derepFs) <- sample.names
names(derepRs) <- sample.names

dadaFs <- dada(derepFs, err=errF, multithread=TRUE, pool="pseudo")
dadaRs <- dada(derepRs, err=errR, multithread=TRUE, pool="pseudo")
# 检查，这里是检测正向的第一个样本
dadaFs[[1]]

mergers <- mergePairs(dadaFs, derepFs, dadaRs, derepRs, trimOverhang=TRUE, minOverlap=170)

seqtab <- makeSequenceTable(mergers)
dim(seqtab)
# 检测长度分布
table(nchar(getSequences(seqtab)))

seqtab.nochim <- removeBimeraDenovo(seqtab, method="consensus", multithread=T, verbose=T)
dim(seqtab.nochim)

sum(seqtab.nochim)/sum(seqtab)

taxa <- assignTaxonomy(seqtab.nochim, "Database/silva_nr_v132_train_set.fa.gz", multithread=T, tryRC=T)

## making and writing out standard output files:
# giving our seq headers more manageable names (ASV_1, ASV_2...)
asv_seqs <- colnames(seqtab.nochim)
asv_headers <- vector(dim(seqtab.nochim)[2], mode="character")

for (i in 1:dim(seqtab.nochim)[2]) {
  asv_headers[i] <- paste(">ASV", i, sep="_")
}

# fasta:
asv_fasta <- c(rbind(asv_headers, asv_seqs))
write(asv_fasta, "DADA2/ASVs.fa")

# count table:
asv_tab <- t(seqtab.nochim)
row.names(asv_tab) <- sub(">", "", asv_headers)
write.table(asv_tab, "DADA2/ASVs_counts.txt", sep="\t", quote=F)

# tax table:
asv_tax <- taxa
row.names(asv_tax) <- sub(">", "", asv_headers)
write.table(asv_tax, "DADA2/ASVs_taxonomy.txt", sep="\t", quote=F)
