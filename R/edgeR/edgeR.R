library(edgeR)
library(statmod)

sampleNames <- c("siCtrl_1", "siCtrl_2", "siSUZ12_1", "siSUZ12_2")
# 第一行是命令信息，所以跳过
data <- read.table("all_feature.txt", header=TRUE, quote="\t", skip=1)
names(data)[7:10] <- sampleNames
countData <- as.matrix(data[7:10])
rownames(countData) <- data$Geneid
condition <- factor(c("siCtrl", "siCtrl", "siSUZ12", "siSUZ12"))

# make edgeR object
genelist <- DGEList(counts=countData, group=condition)

# filter set 0.4 = min(readcounts)/1e6
keep <- rowSums(cpm(genelist) > 0.4) >= 2
table(keep)
genelist.filted <- genelist[keep, , keep.lib.sizes=FALSE]

#
genelist.norm <- calcNormFactors(genelist.filted)

# design matrix
design <- model.matrix(~0+condition)
colnames(design) <- levels(condition)
design

# desprison
genelist.disp <- estimateDisp(genelist.norm, design, robust=TRUE)
plotBCV(genelist.disp)

# quasi-likelihood
fit <- glmQLFit(genelist.disp, design, robust=TRUE)

###
edg <- makeContrasts(siCtrl-siSUZ12, levels=design)
res <- glmQLFTest(fit, contrast=edg)
ig.edger <- res$table[p.adjust(res$table$PValue, method="BH") < 0.1, ]
res$table$padj <- p.adjust(res$table$PValue, method="BH")
write.csv(res$table, "edg_output.csv")
write.csv(merge(res$table, countData, by="row.names", sort=FALSE), "all_edg_output.csv", row.names=FALSE)


# plot
topTags(res, n=10)
is.de <- decideTestsDGE(res)
summary(is.de)
plotMD(res, status=is.de, values=c(1, -1), col=c("red", "green"), legend="topright")

tr <- glmTreat(fit, contrast=edg, lfc=log2(1.5))













