# pzw
# 20180920

# 载入qqman包
library("qqman")

# 读入
SCZ <- read.table('SNP/rall.txt', header=TRUE, sep="\t")

# 查看数据结构
str(SCZ)

# 查看曼哈顿图的帮助说明
?manhattan

# 查看hg19chrc这一列的值（这里会自动去重）
levels(SCZ$hg19chrc)

# 由于曼哈顿图需要的CHR变量是数字(int)，所以需要先把原始文件里的hg19chrc这一列改造一下
# 如把chr1改成1，把chrX改成23
# 去掉chr
SCZ$hg19chrc <- gsub("chr", "", SCZ$hg19chrc)

# 简单的看看去除的效果，看头尾
head(SCZ)
tail(SCZ)

# 将X改成23
SCZ$hg19chrc <- gsub("X", "23", SCZ$hg19chrc)
tail(SCZ)

# 查看当前类型是否int
is.integer(SCZ$hg19chrc)

# 把hg19chrc这一列的类型改成int
SCZ$hg19chrc <- as.integer(SCZ$hg19chrc)

# 再看一次
is.integer(SCZ$hg19chrc)

# 修改列名
colnames(SCZ)[1] <- "CHR"
colnames(SCZ)[5] <- "BP"
colnames(SCZ)[9] <- "P"

# 查看数据结构
str(SCZ)

# 大概是需要一个基准？
# 其实可以这样做，然后下面把genomeWideLine填进-log10里面，这样更优雅
genomeWideLine = 0.05 / nrow(SCZ)

# 画曼哈顿图
manhattan(SCZ, chrlabs=c(1:22, "X"), suggestiveline=FALSE,
          genomewideline=-log10(5.3e-09), main="PGC SCZ GWAS")

# 说要输出最明显的那个染色体，明显是6号咯
manhattan(subset(SCZ, CHR==6), genomewideline=-log10(5.3e-09), snp="snpid",
          suggestiveline=FALSE, main="PGC SCZ GWAS chr6")

# 然后把pvalue最小的几个弄出来，哪个最显著也是很明显的
# rs115329265
manhattan(subset(SCZ, CHR==6), genomewideline=-log10(5.3e-09), suggestiveline=FALSE,
          snp = "snpid",
          main="significant locus", annotatePval=1e-30, annotateTop=TRUE)



