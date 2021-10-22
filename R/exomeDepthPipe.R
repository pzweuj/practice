#!/usr/bin/env Rscript
# 20211012
# pzw

library(ExomeDepth)
library(getopt)
# data(Conrad.hg19)
# data(exons.hg19)

# 传参
# args <- commandArgs(TRUE)
spec <- matrix(
    c(
        "bed", "b", 1, "character", "bed file, need 4 columns without header",
        "input", "i", 1, "character", "bam file directory",
        "output", "o", 1, "character", "output directory"
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

# 输出文件夹
outputPath <- opt$output
if(!(dir.exists(outputPath))) {
    dir.create(outputPath)
}

# 获得counts
my.counts <- getBamCounts(bed.frame=bedFile,
    bam.files=bamFile,
    include.chr=FALSE
)
my.counts.dafr <- as(my.counts, "data.frame")
exomeCount.mat <- as.matrix(my.counts.dafr[, grep(names(my.counts.dafr), pattern="*.bam")])

# 注释
# genes.GRanges.hg19 <- GenomicRanges::GRanges(
#     seqnames=exons.hg19$chromosome,
#     IRanges::IRanges(start=exons.hg19$start, end=exons.hg19$end),
#     names=exons.hg19$name
# )

# 循环运行
nSamples <- ncol(exomeCount.mat)
for (i in 1:nSamples) {
    my.choice <- select.reference.set(
        test.counts=exomeCount.mat[, i],
        reference.counts=exomeCount.mat[, -i],
        bin.length=(my.counts.dafr$end - my.counts.dafr$start)/1000, n.bins.reduced=10000
    )
    my.reference.selected <- apply(X=exomeCount.mat[, my.choice$reference.choice, drop=FALSE], MAR=1, FUN=sum)
    all.exons <- new("ExomeDepth", test=exomeCount.mat[, i],
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
    # all.exons <- AnnotateExtra(x=all.exons, reference.annotation=Conrad.hg19.common.CNVs,
    #     min.overlap=0.5, column.name="Conrad.hg19"
    # )
    # all.exons <- AnnotateExtra(x=all.exons, reference.annotation=genes.GRanges.hg19,
    #     min.overlap=0.0001, column.name="exons.hg19"
    # )
    output.file <- paste(outputPath, "/", colnames(exomeCount.mat)[i], ".txt", sep="")
    write.table(file=output.file, x=all.exons@CNV.calls, row.names=FALSE, quote=FALSE, sep="\t")
}
