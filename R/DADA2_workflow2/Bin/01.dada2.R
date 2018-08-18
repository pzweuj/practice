library(dada2) # loads DADA2

list.files("Rawdata")

# setting a few variables we're going to use
fnFs <- sort(list.files("Rawdata", pattern="_R1_001.fastq.gz", full.names=TRUE))
fnRs <- sort(list.files("Rawdata", pattern="_R2_001.fastq.gz", full.names=TRUE))

sample.names <- sapply(strsplit(basename(fnFs), "_"), `[`, 1)

filtFs <- file.path("Cleandata", paste0(sample.names, "_R1_filtered.fq.gz"))
filtRs <- file.path("Cleandata", paste0(sample.names, "_R2_filtered.fq.gz"))

plotQualityProfile(fnFs[1])
plotQualityProfile(fnRs)
plotQualityProfile(fnFs[17:20])


## 过滤
out <- filterAndTrim(
  fnFs, filtFs, fnRs, filtRs, truncLen=c(240, 160),
  maxN=0, maxEE=c(2, 2), truncQ=2, rm.phix=TRUE, 
  compress=TRUE, multithread=FALSE # 在windows下，multithread设置成FALSE
)
plotQualityProfile(filtFs[1])

## 计算错误率模型
# 分别计算正向和反向
errF <- learnErrors(filtFs, multithread=TRUE)
errR <- learnErrors(filtRs, multithread=TRUE)
# 画出错误率统计图
plotErrors(errF, nominalQ=TRUE)
plotErrors(errR, nominalQ=TRUE)

##消除误差
# 下面其实是一个批量的操作，如果是处理大文件，内存可能不足，更好的做法是一个一个样本的进行
derepFs <- derepFastq(filtFs, verbose=TRUE)
derepRs <- derepFastq(filtRs, verbose=TRUE)
# 用sample.names来改名
names(derepFs) <- sample.names
names(derepRs) <- sample.names
## dada2核心算法
# 从头OTU方法必须在处理之前对样本进行聚类，因为样本之间没有聚类标签不一致且无法比较，即样本1中的OTU1和样本2中的OTU1可能不相同。
# 而DADA2可以更精确地解析序列变异，可以独立处理然后组合。
dadaFs <- dada(derepFs, err=errF, multithread=TRUE)
dadaRs <- dada(derepRs, err=errR, multithread=TRUE)

## 合并双端
mergers <- mergePairs(dadaFs, derepFs, dadaRs, derepRs, verbose=TRUE)
# 构造列表
seqtab <- makeSequenceTable(mergers)
dim(seqtab)
# 检测长度分布
table(nchar(getSequences(seqtab)))

## 去除嵌合体
seqtab.nochim <- removeBimeraDenovo(seqtab, method="consensus", multithread=T, verbose=T)
dim(seqtab.nochim)
sum(seqtab.nochim)/sum(seqtab)

## 物种分类
# 这里是使用silva数据库进行注释。
# 然后把特征序列都提出来，在对这些序列进行命名，命名的方式为“ASV_x”。
taxa <- assignTaxonomy(seqtab.nochim, "Database/silva_nr_v132_train_set.fa.gz", multithread=T, tryRC=T)

asv_seqs <- colnames(seqtab.nochim)
asv_headers <- vector(dim(seqtab.nochim)[2], mode="character")

for (i in 1:dim(seqtab.nochim)[2]) {
  asv_headers[i] <- paste(">ASV", i, sep="_")
}

## 输出DADA2结果
# 这一步输出三个文件，一个是把特征序列都放在一起的fasta文件，一个是注释文件，一个是counts数文件。
# fasta:
asv_fasta <- c(rbind(asv_headers, asv_seqs))
write(asv_fasta, "ASV/ASVs.fa")

# count table:
asv_count <- t(seqtab.nochim)
row.names(asv_count) <- sub(">", "", asv_headers)
write.table(asv_count, "ASV/ASVs_counts.txt", sep="\t", quote=F)

# tax table:
asv_taxa <- taxa
row.names(asv_taxa) <- sub(">", "", asv_headers)
write.table(asv_taxa, "ASV/ASVs_taxonomy.txt", sep="\t", quote=F)

quit("no")
