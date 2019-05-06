# coding=utf-8
# pzw
# 20190505
# v0.1

import sys
import argparse
import pandas


# 得到升高或降低的读数
def getUpLocationID(dataframe, sample, maxCutOff):
	upList = list(dataframe[dataframe[sample] > maxCutOff].index)
	return upList

def getDownLocationID(dataframe, sample, minCutOff):
	downList = list(dataframe[dataframe[sample] < minCutOff].index)
	return downList


def main(inputFile, outputFile, bedFile, minCutOff, maxCutOff, del_length):

	# 读入划分的bed文件
	# sampleID = inputFile.split(".txt")[0]
	# if sampleID.__contains__("/"):
	# 	sampleID = sampleID.split("/")[1]
	bedFile = pandas.read_csv(bedFile, header=None, sep="\t")
	bedFile.rename(columns={0:"chromosome", 1:"start", 2:"end", 3:"gene", 4:"ID", 5:"length"}, inplace=True)

	# 读入上一步的输出结果
	df = pandas.read_csv(inputFile, header=0, sep="\t")
	df.rename(columns={"Unnamed: 0":"ID"}, inplace=True)
	sampleID = df.columns.values.tolist()[1]

	# 获得升高或降低的dataframe
	upDF = bedFile.iloc[getUpLocationID(df, sampleID, maxCutOff), :]
	downDF = bedFile.iloc[getDownLocationID(df, sampleID, minCutOff), :]

	# 处理dataframe，合并相邻区域，这里存在warning，未找到解决方法，不影响结果
	i = 0
	indexList = []
	while i < len(downDF.index) - 1:
		if downDF.index[i+1] - downDF.index[i] == 1:
			start_tmp = downDF.at[downDF.index[i], "start"]
			downDF.loc[downDF.index[i+1], "start"] = start_tmp
		else:
			indexList.append(downDF.index[i])
		i += 1
	indexList.append(downDF.index[-1])

	# 提取最终的结果
	finalDownDF = downDF.loc[indexList, :]
	for m in finalDownDF.index:
		if finalDownDF.at[m, "end"] - finalDownDF.at[m, "start"] < del_length:
			finalDownDF.drop(m, inplace=True)

	# 输出结果
	finalDownDF[["chromosome", "start", "end", "gene"]].to_csv(outputFile, sep="\t", index=False)

	print("task done")

if __name__ == "__main__":
	parser = argparse.ArgumentParser(description="get long del results",
		prog="SV_Long_Del_v0.1.py",
		usage="python3 SV_Long_Del_v0.1.py -i <input file> -o <output file> -b <bedfile> -min <min cut off> -max <max cut off> -l <del length>")
	parser.add_argument("-v", "--version", action="version", version="Version 0.1 20190505")
	parser.add_argument("-i", "--input", type=str, help="Input the file from SV analysis e.g: XXX.txt")
	parser.add_argument("-o", "--output", type=str, help="Output the result file e.g: XXX.results.txt")
	parser.add_argument("-b", "--bedfile", type=str, help="Input the bedfile which was used to get counts e.g: SV_RC.bed")
	parser.add_argument("-min", "--mincutoff", type=float, help="the min number of the refer region")
	parser.add_argument("-max", "--maxcutoff", type=float, help="the max number of the refer region")
	parser.add_argument("-l", "--length", type=int, help="the min length of the del region")
	if len(sys.argv[1:]) == 0:
		parser.print_help()
		parser.exit()
	args = parser.parse_args()
	main(inputFile=args.input, outputFile=args.output, bedFile=args.bedfile, minCutOff=args.mincutoff, maxCutOff=args.maxcutoff, del_length=args.length)

