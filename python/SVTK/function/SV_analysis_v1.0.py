# coding=utf-8
# pzw
# 20190430
# v1.0

import sys
import argparse
import pandas
import numpy

# UR值为【单样本单位置reads数】除以【单样本总reads数】
def getUR(DF, ID, sample):
	mapped_ID = DF.loc[ID, sample]
	mapped_all = DF[sample].sum()
	UR_sample_ID = float(mapped_ID) / float(mapped_all)
	return UR_sample_ID

# UR_mean值为【多样本同一位置的UR平均值】
def getURmean(DF, ID):
	sampleName = DF.columns.values.tolist()
	sampleName.pop(0)
	UR_all = 0
	for i in sampleName:
		UR_sample_ID = getUR(DF, ID, i)
		UR_all += UR_sample_ID
	UR_mean_ID = UR_all / len(sampleName)
	return UR_mean_ID

# SD值为【多样本同一位置的UR值标准差】
def getSD(DF, ID):
	sampleName = DF.columns.values.tolist()
	sampleName.pop(0)
	sqr_sum = 0
	for i in sampleName:
		UR_sample_ID = getUR(DF, ID, i)
		UR_gap_sample_ID = UR_sample_ID - getURmean(DF, ID)
		sqr_sum += (UR_gap_sample_ID ** 2)
	SD_ID = numpy.sqrt(sqr_sum / len(sampleName))
	return SD_ID

# Z值为【单样本单个位置作Z检验得到的Z值】
def getZscore(DF, ID, sample):
	UR_sample_ID = getUR(DF, ID, sample)
	UR_mean_ID = getURmean(DF, ID)
	SD_ID = getSD(DF, ID)
	Z_sample_ID = (UR_sample_ID - UR_mean_ID) / SD_ID
	return Z_sample_ID

# DQ值为【单样本单个位置的UR值】比上【多样本相同位置的UR值的平均值】
def getDQ(DF, ID, sample):
	DQ = getUR(DF, ID, sample) / getURmean(DF, ID)
	return DQ

# 对测试样本标准化
def testNormalize(sampleName):
	sampleIDs = sampleName.split("/")
	for i in sampleIDs:
		if i.__contains__(".counts"):
			sampleID = i.split(".")[0]
	test = open(sampleName, "r")
	n = 1
	ID_List = []
	counts_List = []
	for line in test:
		if line.startswith("#"):
			continue
		else:
			lineAS = line.split("\t")
			ID_List.append("ID_" + str(n).zfill(3))
			counts_List.append(float(lineAS[3]))
			n += 1
	sample_counts = {"ID": ID_List, sampleID: counts_List}
	testSample = pandas.DataFrame(sample_counts)
	testSample = testSample[["ID", sampleID]]
	test.close()
	return testSample

def main(referenceFile, testSample, testMethod, outputFile):
	reference = pandas.read_csv(referenceFile, header=0, sep="\t")
	reference.rename(columns={"Unnamed: 0":"ID"}, inplace=True)
	df = testNormalize(testSample).merge(reference, how="inner", on="ID")
	df.index = df["ID"].tolist()
	del reference
	del df["ID"]
	
	# 创建新表
	colNamesList = df.columns.values.tolist()
	indexNamesList = df.index.tolist()
	df_new = pandas.DataFrame(index=indexNamesList, columns=colNamesList)
	
	cal = 0
	for colNames in colNamesList:
		for indexNames in indexNamesList:
			cal += 1
			cellAmount = df.shape[0] * df.shape[1]
			print("process: " + "%.2f" % ((float(cal) / cellAmount) * 100) + "%", end="\r")
			if testMethod == "Z":
				df_new.at[indexNames, colNames] = "%.2f" % getZscore(df, indexNames, colNames)
			elif testMethod == "DQ":
				df_new.at[indexNames, colNames] = "%.2f" % getDQ(df, indexNames, colNames)
			else:
				print("please input 'Z' or 'DQ' to chose analysis method!")

	df_new.to_csv(outputFile, sep="\t")
	print("\n")
	print("task done")

if __name__ == "__main__":
	parser = argparse.ArgumentParser(description="calculate Z score or DQ score",
		prog="SV_analysis_v1.0.py",
		usage="python3 SV_analysis_v1.0.py -r <reference file> -i <test file> -m <method> -o <output file>")
	parser.add_argument("-v", "--version", action="version", version="Version 1.0 20190430")
	parser.add_argument("-r", "--reference", type=str, help="Input the the reference files e.g: reference.txt")
	parser.add_argument("-i", "--input", type=str, help="Input the the test file e.g: XXX.counts")
	parser.add_argument("-m", "--method", type=str, help="chose a method e.g: Z or DQ")
	parser.add_argument("-o", "--output", type=str, help="output the results e.g: XXX.DQ.txt")
	if len(sys.argv[1:]) == 0:
		parser.print_help()
		parser.exit()
	args = parser.parse_args()
	main(referenceFile=args.reference, testSample=args.input, testMethod=args.method, outputFile=args.output)
