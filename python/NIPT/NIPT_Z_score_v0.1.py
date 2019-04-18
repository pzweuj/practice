# encoding:utf-8
# pzw
# 20190320
# v0.1 已完成Z检验
###########
import pandas
import numpy
import sys

##############
testSample = pandas.read_csv(sys.argv[1], names=["chrom", "test"], sep="\t")
reference = pandas.read_csv(sys.argv[1], header=0, sep="\t")
##############

### test #####
# testSample = pandas.read_csv("H05078S.counts", names=["chrom", "test"], sep="\t")
# reference = pandas.read_csv("reference_counts.txt", header=0, sep="\t")
##############

###########
# 预处理 将受检数据合并到参考数据
reference.rename(columns={"Unnamed: 0":"chrom"}, inplace=True)
df = reference.merge(testSample, how="inner", on="chrom")
del testSample
del reference
df.index = df["chrom"].tolist()
del df["chrom"]

# 定义UR值函数
def getUR(chrom, sample):
	mapped_chrom = df.loc[chrom, sample]
	mapped_all = df[sample].sum()
	UR_sample_chrom = float(mapped_chrom) / float(mapped_all)
	return UR_sample_chrom

# 定义UR均值函数
def getURmean(chrom):
	sampleName = df.columns.values.tolist()
	sampleName.pop(0)
	UR_all = 0
	for i in sampleName:
		UR_sample_chrom = getUR(chrom, i)
		UR_all += UR_sample_chrom
	UR_mean_chrom = UR_all / len(sampleName)
	return UR_mean_chrom

# 定义标准差函数
def getSD(chrom):
	sampleName = df.columns.values.tolist()
	sampleName.pop(0)
	sqr_sum = 0
	for i in sampleName:
		UR_sample_chrom = getUR(chrom, i)
		UR_gap_sample_chrom = UR_sample_chrom - getURmean(chrom)
		sqr_sum += UR_gap_sample_chrom ** 2
	SD_chrom = numpy.sqrt(sqr_sum / len(sampleName))
	return SD_chrom

# 定义Z检验函数
def getZscore(chrom, sample):
	UR_sample_chrom = getUR(chrom, sample)
	UR_mean_chrom = getURmean(chrom)
	SD_chrom = getSD(chrom)
	Z_sample_chrom = (UR_sample_chrom - UR_mean_chrom) / SD_chrom
	return Z_sample_chrom

# 获得Z值
risk_13 = "%.2f" % getZscore("chr13", "test")
risk_18 = "%.2f" % getZscore("chr18", "test")
risk_21 = "%.2f" % getZscore("chr21", "test")

print "13号染色体三体风险值为： ", risk_13
print "18号染色体三体风险值为： ", risk_18
print "21号染色体三体风险值为： ", risk_21