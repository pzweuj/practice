library(WGCNA)

# 设置将strings不要转成factors
options(stringsAsFactors=F)

# 读入表达数据
femData <- read.csv("LiverFemale3600.csv")
# 查看数据
dim(femData)
names(femData)

# 数据初步处理
# 提取出表达数据
datExpr0 <- as.data.frame(t(femData[, -c(1:8)]))
names(datExpr0) <- femData$substanceBXH
rownames(datExpr0) <- names(femData)[-c(1:8)]

# 样本聚类检查离群值（就是树上特别不接近的）
gsg <- goodSamplesGenes(datExpr0, verbose=3)
# 是true的话说明所有genes都通过了筛选
gsg$allOK

if (!gsg$allOK){
  if (sum(!gsg$goodGenes)>0)
    printFlush(paste("Removing genes:", paste(names(datExpr0)[!gsg$goodGenes], collapse=", ")))
  if (sum(!gsg$goodSamples)>0) 
    printFlush(paste("Removing samples:", paste(rownames(datExpr0)[!gsg$goodSamples], collapse=", ")))
  datExpr0 <- datExpr0[gsg$goodSamples, gsg$goodGenes]
}

# 画图来看是否有特别离群的
sampleTree <- hclust(dist(datExpr0), method="average")
sizeGrWindow(12, 9)
par(cex=0.6)
par(mar=c(0,4,2,0))
plot(sampleTree, main="Sample clustering to detect outliers", sub="", xlab="", cex.lab=1.5, 
     cex.axis=1.5, cex.main=2)
abline(h=15, col="red")

# 把离群的修剪掉
clust <- cutreeStatic(sampleTree, cutHeight=15, minSize=10)
table(clust)
keepSamples <- (clust==1)
datExpr <- datExpr0[keepSamples, ]
nGenes <- ncol(datExpr)
nSamples <- nrow(datExpr)

# 读入性状数据
traitData <- read.csv("ClinicalTraits.csv")
dim(traitData)
names(traitData)

# 只保留需要的信息
allTraits <- traitData[, -c(31, 16)]
allTraits <- allTraits[, c(2, 11:36)]
dim(allTraits)
names(allTraits)

# 使两个数据框的行名一致
femaleSamples <- rownames(datExpr)
traitRows <- match(femaleSamples, allTraits$Mice)
datTraits <- allTraits[traitRows, -1]
rownames(datTraits) <- allTraits[traitRows, 1]

# 清理，释放内存
collectGarbage()

# 作性状与样本的关系热图
sampleTree2 <- hclust(dist(datExpr), method="average")
traitColors <- numbers2colors(datTraits, signed=F)
plotDendroAndColors(sampleTree2, traitColors,
                    groupLabels=names(datTraits), 
                    main="Sample dendrogram and trait heatmap")

# 保存预处理完成的数据
save(datExpr, datTraits, file="dataInput.RData")

# 允许使用最大线程
allowWGCNAThreads()
# 或者直接指定线程数
enableWGCNAThreads(nThreads=16)

# 软阈值的预设范围
powers <- c(c(1:10), seq(from=12, to=20, by=2))
# 自动计算推荐的软阈值
sft <- pickSoftThreshold(datExpr, powerVector=powers, verbose=5, networkType="unsigned")
# 推荐值。如果是NA，就需要画图来自己挑选
sft$powerEstimate

# 作图
sizeGrWindow(9, 5)
par(mfrow=c(1, 2))
cex1 <- 0.9
plot(sft$fitIndices[,1], -sign(sft$fitIndices[, 3])*sft$fitIndices[, 2],
     xlab="Soft Threshold (power)", ylab="Scale Free Topology Model Fit,signed R^2", type="n",
     main=paste("Scale independence"))
text(sft$fitIndices[, 1], -sign(sft$fitIndices[, 3])*sft$fitIndices[, 2],
     labels=powers, cex=cex1, col="red")
abline(h=0.80, col="red")
plot(sft$fitIndices[, 1], sft$fitIndices[, 5],
     xlab="Soft Threshold (power)", ylab="Mean Connectivity", type="n",
     main=paste("Mean connectivity"))
text(sft$fitIndices[, 1], sft$fitIndices[, 5], labels=powers, cex=cex1, col="red")

# 构建网络
# deepSplit是生成模块的梯度，从[0:4]中选择，越大模块越多
# minModuleSize最小模块的基因数目
# mergeCutHeight是合并相似的模块的合并系数，是通过主成分分析出来的
# mumericLabels 已数字命名模块
# nThreads 线程数，当设置为0时使用最大线程
net <- blockwiseModules(datExpr, corType="pearson",
                        power=sft$powerEstimate,
                        TOMType="unsigned", saveTOMs=TRUE, saveTOMFileBase="femaleMouseTOM",
                        deepSplit=2, minModuleSize=30,
                        reassignThreshold=0, mergeCutHeight=0.25,
                        numericLabels=T, pamRespectsDendro=F, nThreads=0,
                        verbose=3)

# 查看每个模块的基因数，其中0模块下为没有计算进入模块的基因数
table(net$colors)

# 可视化
sizeGrWindow(12, 9)
# 把模组编号转成颜色
mergedColors <- labels2colors(net$colors)
plotDendroAndColors(net$dendrograms[[1]], mergedColors[net$blockGenes[[1]]],
                    "Module colors",
                    dendroLabels=F, hang=0.03,
                    addGuide=T, guideHang=0.05)
# 计算模块特征向量MEs
moduleLabels <- net$colors
moduleColors <- labels2colors(net$colors)
MEs <- net$MEs
geneTree <- net$dendrograms[[1]]
# 保存数据
save(MEs, moduleLabels, moduleColors, geneTree,
     file="networkConstruction-auto.RData")


nGenes <- ncol(datExpr)
nSamples <- nrow(datExpr)

# 重新计算MEs
MEs0 <- moduleEigengenes(datExpr, moduleColors)$eigengenes
MEs <- orderMEs(MEs0)
# 主成分向量和性状的关联，pearson校正
moduleTraitCor <- cor(MEs, datTraits, use="p")
moduleTraitPvalue <- corPvalueStudent(moduleTraitCor, nSamples)

# 画出模组和性状的关联热图
sizeGrWindow(10,6)
# 把校正值和p值写在一起
textMatrix <-  paste(signif(moduleTraitCor, 2), "\n(",
                     signif(moduleTraitPvalue, 1), ")", sep="")
dim(textMatrix) = dim(moduleTraitCor)
par(mar=c(6, 8.5, 3, 3))
# 画图
labeledHeatmap(Matrix=moduleTraitCor,
               xLabels=names(datTraits),
               yLabels=names(MEs),
               ySymbols=names(MEs),
               colorLabels=F,
               colors=greenWhiteRed(50),
               textMatrix=textMatrix,
               setStdMargins=F,
               cex.text=0.5,
               zlim=c(-1,1),
               main=paste("Module-trait relationships"))

## 计算模组
modNames <- substring(names(MEs), 3)
# 计算MM
geneModuleMembership <- as.data.frame(cor(datExpr, MEs, use="p"))
MMPvalue <- as.data.frame(corPvalueStudent(as.matrix(geneModuleMembership), nSamples))
names(geneModuleMembership) <- paste("MM", modNames, sep="")
names(MMPvalue) <- paste("p.MM", modNames, sep="")
# 计算GS
geneTraitSignificance <- as.data.frame(cor(datExpr, datTraits, use="p"))
GSPvalue <- as.data.frame(corPvalueStudent(as.matrix(geneTraitSignificance), nSamples))
names(geneTraitSignificance) <- paste("GS.", names(datTraits), sep="")
names(GSPvalue) <- paste("p.GS.", names(datTraits), sep="")

# 输出
geneInfo <- cbind(geneModuleMembership, MMPvalue, geneTraitSignificance, GSPvalue)
write.table(geneInfo, file="geneInfo.txt", sep="\t", quote=F)


# power是之前的软阈值
TOM <- TOMsimilarityFromExpr(datExpr, power=sft$powerEstimate)
dissTOM <- 1 - TOM
#为了更显著，用7次方
plotTOM <- dissTOM^7
diag(plotTOM) <- NA
TOMplot(plotTOM, geneTree, moduleColors, main="Network heatmap plot, all genes")


# 随机抽，自行取值
nSelect <- 400
set.seed(10)
select <- sample(nGenes, size=nSelect)
selectTOM <- dissTOM[select, select]

# 重新聚类
selectTree <- hclust(as.dist(selectTOM), method="average")
selectColors <- moduleColors[select]
# 绘图
plotDiss = selectTOM^7
diag(plotDiss) = NA
TOMplot(plotDiss, selectTree, selectColors, main="Network heatmap plot, selected genes")
# 模块聚类与热图
plotEigengeneNetworks(MEs, "", marDendro=c(0,4,1,2),
                      marHeatmap=c(3,4,1,2), cex.lab=0.8,
                      plotDendrograms=TRUE, plotHeatmaps=TRUE,
                      xLabelsAngle=90)


# cytoscape
# 选择可视化模块
modules <- c("brown")
# 得到基因ID
probes <- colnames(datExpr)
inModule <- is.finite(match(moduleColors, modules))
modProbes <- probes[inModule]
# 得到候选基因的TOM
modTOM <- TOM[inModule, inModule]
# 命名
dimnames(modTOM) <- list(modProbes, modProbes)
# 输出cytoscape
cyt <- exportNetworkToCytoscape(modTOM,
                                edgeFile=paste("brown-cyto-edge-", paste(modules, collapse="-"), ".txt", sep=""),
                                nodeFile=paste("brown-cyto-nodes-", paste(modules, collapse="-"), ".txt", sep=""),
                                weighted=T, threshold = 0.02,
                                nodeNames=modProbes,
                                nodeAttr=moduleColors[inModule])




