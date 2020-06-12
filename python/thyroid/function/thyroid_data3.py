# coding=utf-8
# pzw
# 20200521
# 不进行丰度过滤

import os
import sys
import argparse

def getAbsPath():
	now = os.path.abspath(os.path.dirname(sys.argv[0]))
	return now

def GetPanelResults(annovarFileName, bed, depth):
	dbCheck = open(bed, "r")
	results = {}
	for b in dbCheck:
		if b.startswith("chrom"):
			continue
		else:
			bSplit = b.split("\n")[0].split("\t")
			chrom = bSplit[0]
			start = bSplit[1]
			end = bSplit[2]
			uniqueID = bSplit[3]

			results[uniqueID] = "0"
			annovarFile = open(annovarFileName, "r")
			print(annovarFileName)
			for a in annovarFile:
				if a.startswith("Chr\tStart"):
					continue
				else:
					aSplit = a.split("\t")
					chromA = aSplit[0]
					startA = aSplit[1]
					endA = aSplit[2]

					otherinfo = aSplit[19]
					DP = int(otherinfo.split(";")[0].split("=")[1])
					AF = aSplit[21].split(":")[2]
					if "," in AF:
						AF = float(AF.split(",")[0])
					else:
						AF = float(AF)

					if chromA == chrom:
						if startA == start:
							if endA == end:
								if DP >= depth:
									results[uniqueID] = str(AF)
			annovarFile.close()
	dbCheck.close()
	resultsList = []
	for i in sorted(results.keys(), key=int):
		resultsList.append(results[i])
	return resultsList





def AnalysisPath(AnnovarFilePath, outputArff, bed, depth, trainDict=""):
	file_list = os.listdir(AnnovarFilePath)
	test = open(outputArff, "w")

	dbCheck = open(bed, "r")
	n = 0
	for l in dbCheck:
		if l.startswith("chrom"):
			continue
		else:
			n += 1

	i = 0
	numList = []
	while i <= (n - 1):
		numList.append("T" + str(i))
		i += 1
	test.write("SampleName\tClass\t" + "\t".join(numList) + "\n")

	sampleList = []
	for f in file_list:
		resultsList = GetPanelResults(AnnovarFilePath + "/" + f, bed, depth)
		sampleName = f.split(".")[0]
		sampleList.append(sampleName)
		oo = "\t".join(resultsList) + "\n"
		if trainDict:
			test.write(sampleName + "\t" + trainDict[sampleName] + "\t" + oo)
		else:
			test.write(sampleName + "\t" + "?\t" + oo)

	test.close()
	return sampleList

def main(AnnovarFilePath, outputArff, bed, depth, train):
	if train:
		trainDict = {}
		trainFile = open(train, "r")
		for line in trainFile:
			lines = line.split("\t")
			sampleName = lines[0]
			diagnoise = lines[1].split("\n")[0]
			trainDict[sampleName] = diagnoise

		sampleList = AnalysisPath(AnnovarFilePath, outputArff, bed, depth, trainDict)
	else:
		sampleList = AnalysisPath(AnnovarFilePath, outputArff, bed, depth)

	print("Task done")

if __name__ == "__main__":
	parser = argparse.ArgumentParser(description="Arff File maker",
		prog="thyroid_data.py",
		usage="python3 thyroid_data.py -i <input dir> -o <output arff> -b <bed file>")
	group = parser.add_mutually_exclusive_group()
	parser.add_argument("-v", "--version", action="version",
		version="Version 0.1 20200331")
	parser.add_argument("-i", "--input", type=str,
		help="annovar file directory")
	parser.add_argument("-o", "--output", type=str,
		help="the output arff file")
	parser.add_argument("-b", "--bed", type=str,
		help="panel bed file，格式参考function文件夹内bed文件。")
	parser.add_argument("-d", "--depth", type=int,
		help="位点过滤深度，默认500", default=500)
	parser.add_argument("-t", "--train", type=str, default="",
		help="train infomation, option")
	if len(sys.argv[1:]) == 0:
		parser.print_help()
		parser.exit()
	args = parser.parse_args()
	main(AnnovarFilePath=args.input, outputArff=args.output, bed=args.bed, depth=args.depth, train=args.train)
