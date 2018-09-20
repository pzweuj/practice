#pzw
#20180830

#source("https://bioconductor.org/biocLite.R")
#biocLite("DOSE") 
#biocLite("topGO")
#biocLite("clusterProfiler")
#biocLite("pathview")

library(DOSE)
library(org.Hs.eg.db)
library(clusterProfiler)
library(ggplot2)

data <- read.table("gene.list",header=TRUE)
data$GeneName <- as.character(data$GeneName)
transID = bitr(data$GeneName,
               fromType="SYMBOL",
               toType=c("ENSEMBL", "ENTREZID"),
               OrgDb="org.Hs.eg.db"
)

# CC
CC <- enrichGO(transID$ENTREZID,
               "org.Hs.eg.db",
               keyType="ENTREZID",
               ont="CC",
               pvalueCutoff=0.05,
               pAdjustMethod="BH",
               qvalueCutoff=0.1
)
CC <- setReadable(CC, OrgDb=org.Hs.eg.db)
CC_df <- as.data.frame(CC@result)
write.table(CC_df, file="GO_CC.xls", sep="\t", row.names=F)

# MF
MF <- enrichGO(transID$ENTREZID,
               "org.Hs.eg.db",
               keyType="ENTREZID",
               ont="MF",
               pvalueCutoff=0.05,
               pAdjustMethod="BH",
               qvalueCutoff=0.1
)
MF <- setReadable(MF, OrgDb=org.Hs.eg.db)
MF_df <- as.data.frame(MF@result)
write.table(MF_df, file="GO_MF.xls", sep="\t", row.names=F)

# BP
BP <- enrichGO(transID$ENTREZID,
               "org.Hs.eg.db",
               keyType="ENTREZID",
               ont="BP",
               pvalueCutoff=0.05,
               pAdjustMethod="BH",
               qvalueCutoff=0.1
)
BP <- setReadable(BP, OrgDb=org.Hs.eg.db)
BP_df <- as.data.frame(BP@result)
write.table(BP_df, file="GO_BP.xls", sep="\t", row.names=F)

CC_top5 <- head(CC_df, n=5)
MF_top5 <- head(MF_df, n=5)
BP_top5 <- head(BP_df, n=5)

for(i in c(1,2,3,4,5)){
  CC_top5$Description[i] <- paste(CC_top5$Description[i], "CC")
  MF_top5$Description[i] <- paste(MF_top5$Description[i], "MF")
  BP_top5$Description[i] <- paste(BP_top5$Description[i], "BP")
}

GO_df <- rbind(CC_top5, MF_top5, BP_top5)
GO_df$padj <- 0 - GO_df$p.adjust

p <- ggplot(GO_df, aes(x=GO_df$Description, y=GO_df$GeneRatio, fill=GO_df$padj)) + geom_bar(stat="identity") +
  theme(axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1)) +
  ylab("sample number/background number") +
  xlab("pathway name") +
  labs(fill="p.adjust") +
  scale_fill_gradient(low="#C1FFC1",
                      high="#228B22",
                      labels=c("0.01", "0"),
                      breaks=c(-0.01, 0),
                      limits=c(-0.01, 0)
                      )