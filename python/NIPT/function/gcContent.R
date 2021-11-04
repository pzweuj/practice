#
args <- commandArgs(T)
rawCounts <- args[1]
corCounts <- args[2]

RC_DT <- read.table(rawCounts, sep="\t", head=TRUE)
gcCount.loess <- loess(
    RC~GC,
    data=RC_DT,
    control=loess.control(surface="direct"),
    degree=2
    )
predictions <- predict(gcCount.loess, RC_DT$GC)
# resi <- RC_DT$RC - predictions
# RC_DT$RC <- resi
RC_DT$RC <- predictions

write.table(RC_DT, corCounts, sep="\t", quote=FALSE, row.names=FALSE)

