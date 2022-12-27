#!/usr/bin/env Rscript
# 20221226
# pzw

library(ExomeDepth)
library(getopt)

# 传参
# args <- commandArgs(TRUE)
spec <- matrix(
    c(
        "bed", "b", 1, "character", "bed file, need 4 columns without header",
        "input", "i", 1, "character", "bam file directory",
        "output", "o", 1, "character", "output Rdata file"
    ),
    byrow=TRUE, ncol=5   
)
opt <- getopt(spec)
if (is.null(opt$input) || is.null(opt$output) || is.null(opt$bed)) {
    cat(paste(getopt(spec, usage=TRUE), "\n"))
    q()
}

# 导入bed文件
bedFilePath <- opt$bed
bedFile <- read.table(bedFilePath, head=FALSE)
names(bedFile) <- c("chromosome", "start", "end", "name")
bedFile <- as.data.frame(bedFile)

# 导入bam文件
bamFilePath <- opt$input
bamFile <- list.files(bamFilePath, pattern="*.bam$")
bamFile <- file.path(bamFilePath, bamFile)

# 获得counts
my.counts <- getBamCounts(bed.frame=bedFile,
    bam.files=bamFile,
    include.chr=FALSE
)
my.counts.dafr <- as(my.counts, "data.frame")
exomeCount.mat <- as.matrix(my.counts.dafr[, grep(names(my.counts.dafr), pattern="*.bam")])

# 参考导出
referenceData <- opt$output
outputDir <- dirname(opt$output)
if(!(dir.exists(outputDir))) {
    dir.create(outputDir)
}
save(exomeCount.mat, file=referenceData)

