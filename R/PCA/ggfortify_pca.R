# pzw
# 20171203
# coding:utf-8

setwd('~/r/pca')
library(ggfortify)

# 输入文件，这里加了一个check.names=F，因为R语言table不能用数字开头，会自动给你加X的，所以加上这个才不会出错。
final <- read.table('finaldata.txt', header = TRUE, sep = '\t', check.names = FALSE)

# 这些只要样本名字就好了
asian <- read.table('asian.txt', header = FALSE)
black <- read.table('black.txt', header = FALSE)
american <- read.table('american indian.txt', header = FALSE)
white <- read.table('white.txt', header = FALSE)

# 把gene这一列定义为行名
row.names(final) <- final$gene

# 去掉gene这一列
final <- final[-c(1)]

# 行列转置
final <- t(final)

# 定义为dataframe
final <- as.data.frame(final)

# 把final的行名定义为Samples
final$Samples <- rownames(final)

# 定义一个Group，初始值设为NA
final$Group <- NA

# 如果Samples在asian里就把那一行的Group改成ASIAN，以此类推，都不在的Group改成Unknown
final$Group <- ifelse(final$Samples %in% asian$V1, 'ASIAN',
                      ifelse(final$Samples %in% black$V1, 'BLACK',
                             ifelse(final$Samples %in% american$V1, 'AMERICAN INDIAN',
                                    ifelse(final$Samples %in% white$V1, 'WHITE','UnKnown'))))

# 定义一个final.data
final.data <- final[c(1:60483)]

# 把空值用0代替，这里没用的，因为貌似源文件没有空值
# final.data[is.na(final.data)] <- 0

# 画图
# autoplot(prcomp(final.data))

# 画加颜色的图
autoplot(prcomp(final.data), data = final, colour = 'Group')

# 只画那几个(弄成新的表)
# final2 <- subset(final, final$Group == 'ASIAN' | final$Group == 'BLACK' | final$Group == 'AMERICAN INDIAN' | final$Group == 'WHITE' )
# subset.final.data <- final2[c(1:60483)]
# autoplot(prcomp(final.data), data = final2, colour = 'Group')

PCA = prcomp(final.data)
summary <- summary(PCA)
write.table(summary$importance, file = 'HapMap_PCA_summary.txt', sep = '\t', col.names = NA)