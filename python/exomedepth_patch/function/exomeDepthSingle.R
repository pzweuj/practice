#!/usr/bin/env Rscript
# 20221227
# pzw

library(ExomeDepth)
library(getopt)


# 传参
args <- commandArgs(TRUE)
spec <- matrix(
  c(
    "bed", "b", 1, "character", "bed file, need 4 columns without header",
    "input", "i", 1, "character", "bam file",
    "output", "o", 1, "character", "output txt",
    "reference", "r", 1, "character", "reference rdata"
  ),
  byrow=TRUE, ncol=5
)
opt <- getopt(spec)
if (is.null(opt$input) || is.null(opt$output) || is.null(opt$bed) || is.null(opt$reference)) {
  cat(paste(getopt(spec, usage=TRUE), "\n"))
  q()
}


# 导入reference
load(opt$reference)


## 测试
# load("my.counts.rdata")

# 导入bed文件
bedFilePath <- opt$bed
bedFile <- read.table(bedFilePath, head=FALSE)
names(bedFile) <- c("chromosome", "start", "end", "name")
bedFile <- as.data.frame(bedFile)

# 导入bam文件
bamFile <- opt$input
my.counts <- getBamCounts(bed.frame=bedFile,
                          bam.files=bamFile,
                          include.chr=FALSE)
my.counts.dafr <- as(my.counts, "data.frame")
my.counts.mat <- my.counts.dafr[, dim(my.counts.dafr)[2]]
my.counts.name <- colnames(my.counts.dafr)[5]

# 合并测试样本和对照组
my.counts.combine.ref <- cbind(my.counts.mat, exomeCount.mat)
colnames(my.counts.combine.ref)[1] <- my.counts.name

my.choice <- select.reference.set(
  test.counts=my.counts.combine.ref[, 1],
  reference.counts=my.counts.combine.ref[, -1],
  bin.length=(my.counts.dafr$end - my.counts.dafr$start)/1000, n.bins.reduced=10000
)
my.reference.selected <- apply(X=my.counts.combine.ref[, my.choice$reference.choice, drop=FALSE], MAR=1, FUN=sum)
all.exons <- new("ExomeDepth", test=my.counts.combine.ref[, 1],
                 reference=my.reference.selected,
                 formula="cbind(test, reference) ~ 1"
)
all.exons <- CallCNVs(x=all.exons,
                      transition.probability=10^-4,
                      chromosome=my.counts.dafr$chromosome,
                      start=my.counts.dafr$start,
                      end=my.counts.dafr$end,
                      name=my.counts.dafr$exon
)
output.file <- opt$output
write.table(file=output.file, x=all.exons@CNV.calls, row.names=FALSE, quote=FALSE, sep="\t")

# end
