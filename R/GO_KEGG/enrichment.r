#pzw
#20180503

#source("https://bioconductor.org/biocLite.R")
#biocLite("DOSE") 
#biocLite("topGO")
#biocLite("clusterProfiler")
#biocLite("pathview")

library(DOSE)
library(org.Hs.eg.db)
library(topGO)
library(clusterProfiler)
library(pathview)

# import data
data <- read.table("gene.list",header=TRUE)
data$GeneName <- as.character(data$GeneName)
transID = bitr(data$GeneName,
	fromType="SYMBOL",
	toType=c("ENSEMBL", "ENTREZID"),
	OrgDb="org.Hs.eg.db"
)

dir.create("GO")
dir.create("KEGG")

# GO_CC
CC <- enrichGO(transID$ENTREZID,
	"org.Hs.eg.db",
	keyType="ENTREZID",
	ont="CC",
	pvalueCutoff=0.05,
	pAdjustMethod="BH",
	qvalueCutoff=0.1
)
CC <- setReadable(CC, OrgDb=org.Hs.eg.db)

pdf(file="./GO/GO_CC.pdf", bg="transparent")
dotplot(CC, showCategory=12, colorBy="pvalue", font.size=8, title="GO_CC") # + theme(axis.text.y = element_text(angle = 45))
barplot(CC, showCategory=12, title="GO_CC", font.size=8)
plotGOgraph(CC)
dev.off()

write.table(as.data.frame(CC@result), file="./GO/GO_CC.xls", sep="\t", row.names=F)

# GO_MF
MF <- enrichGO(transID$ENTREZID, "org.Hs.eg.db", keyType="ENTREZID", ont="MF", pvalueCutoff=0.05, pAdjustMethod="BH", qvalueCutoff=0.1)
MF <- setReadable(MF, OrgDb=org.Hs.eg.db)

pdf(file="./GO/GO_MF.pdf", bg="transparent")
dotplot(MF, showCategory=12, colorBy="pvalue", font.size=8, title="GO_MF") # + theme(axis.text.y = element_text(angle = 45))
barplot(MF, showCategory=12, title="GO_MF", font.size=8)
plotGOgraph(MF)
dev.off()

write.table(as.data.frame(MF@result), file="./GO/GO_MF.xls", sep="\t", row.names=F)

# GO_BP
BP <- enrichGO(transID$ENTREZID, "org.Hs.eg.db", keyType="ENTREZID", ont="BP", pvalueCutoff=0.05, pAdjustMethod="BH", qvalueCutoff=0.1)
BP <- setReadable(BP, OrgDb=org.Hs.eg.db)

pdf(file="./GO/GO_BP.pdf", bg="transparent")
dotplot(BP, showCategory=12, colorBy="pvalue", font.size=8, title="GO_BP") # + theme(axis.text.y = element_text(angle = 45))
barplot(BP, showCategory=12, title="GO_BP", font.size=8)
plotGOgraph(BP)
dev.off()

write.table(as.data.frame(BP@result), file="./GO/GO_BP.xls", sep="\t", row.names=F)

# KEGG
kegg <- enrichKEGG(transID$ENTREZID, organism="hsa", pvalueCutoff=0.05, pAdjustMethod="BH", qvalueCutoff=0.1)
kegg <- setReadable(kegg, OrgDb=org.Hs.eg.db, keytype="ENTREZID")

pdf(file="./KEGG/KEGG.pdf", bg="transparent")
dotplot(kegg, showCategory=12, colorBy="pvalue", font.size=8, title="KEGG") # + theme(axis.text.y = element_text(angle = 45))
barplot(kegg, showCategory=12, title="KEGG", font.size=8)
dev.off()

write.table(as.data.frame(kegg@result), file="./KEGG/kegg.xls", sep="\t", row.names=F)

dir.create("./KEGG/MAP")
kegg_df = as.data.frame(kegg)

for(i in kegg_df$ID){
  pathview(gene.data=transID$ENTREZID,
           pathway.id=i,
           species="hsa",
           kegg.native=TRUE,
           kegg.dir="./KEGG/MAP"
  )
}

print("TASK DONE")