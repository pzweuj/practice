# coding=utf-8
# pzw
# 20200429

import os
import sys
import argparse
import shutil

def sort_merge(bedfile):
	# 1 排序合并
	bedtools = "bedtools"
	cmd = """
		sed -i 's/\\r//g' {bedfile}
		{bedtools} sort -i {bedfile} | \\
			{bedtools} merge -i - > temp/temp.bed
	""".format(bedfile=bedfile, bedtools=bedtools)
	os.system(cmd)

# 2 补全长度
def fixLength(length):
	tempfile = open("temp/temp.bed", "r")
	checkfile = open("temp/check.bed", "w")
	for line in tempfile:
		cell = line.split("\n")[0].split("\t")
		chrom = cell[0]
		start = cell[1]
		end = cell[2]
		if (int(end) - int(start) + 1) < length:
			fix = length - (int(end) - int(start) + 1)
			fix_start = int(fix / 2)
			if int(start) <= fix_start:
				fix_start = 0
			fix_end = fix - fix_start
			start = str(int(start) - fix_start)
			end = str(int(end) + fix_end)
		output = chrom + ":" + start + "-" + end
		print(output)
		checkfile.write(output + "\n")
	checkfile.close()
	tempfile.close()

def getFasta(outputFile):
	# 3 获得序列
	hg19 = "hg19.fa"
	samtools = "samtools"
	cmd = """
		cat temp/check.bed | while read line; do \\
			{samtools} faidx {hg19} $line; \\
			done > {outputFile}
	""".format(samtools=samtools, hg19=hg19, outputFile=outputFile)
	os.system(cmd)

def Design(outputFasta, outputFile, length, thread, tiling):
	# 4 探针设计
	if tiling <= 0:
		print("[Error]： 叠瓦数只能设置为正整数")
		exit()
	elif isinstance(tiling, float):
		print("[Error]: 叠瓦数只能设置为正整数")
		exit()
	else:
		tileBases = int(length / tiling)
		if tileBases < 3:
			print("[Waring]：叠瓦数设置过大，重定义为1x")
			tileBases = 3
		if tiling == 1:
			tileBases = 3
	length = str(length)
	tileBases = str(tileBases)
	design = "mrbait.py"
	cmd = """
		python3 {design} -A {outputFile} -b {length} -o {outputFasta} \\
			-T {thread} -s tile={tileBases}
	""".format(design=design, outputFile=outputFile, length=length, outputFasta=outputFasta, thread=thread, tileBases=tileBases)
	os.system(cmd)

def rmDuplicates(probeFasta, outputFasta, outputInfo):
	probeFile = open(probeFasta, "r")
	outputFastaFile = open(outputFasta, "w")
	outputInfoFile = open(outputInfo, "w")
	fasta = {}
	k = ""
	inputNum = 0
	for line in probeFile:
		if line.startswith(">"):
			k = line
			fasta[k] = ""
			inputNum += 1
		else:
			fasta[k] = fasta[k] + line.split("\n")[0]
	probeFile.close()

	valueDict = {}
	for i in list(fasta.keys()):
		v = fasta[i].upper()
		if v in valueDict.keys():
			del fasta[i]
			valueDict[v].append(i.split("_Bait")[0].split(">")[1])
		else:
			valueDict[v] = [i.split("_Bait")[0].split(">")[1]]
	outputNum = len(fasta)
	dupsNum = 0
	for j in valueDict.keys():
		if len(valueDict[j]) > 1:
			dupsNum += len(valueDict[j])
	uniqueNum = dupsNum - (inputNum - outputNum)

	outputInfoFile.write("Unique Probes: " + str(outputNum) + "\n")
	outputInfoFile.write("Original Input Probes: " + str(inputNum) + "\n\n")
	outputInfoFile.write("Total Probes Condensed: " + str(dupsNum) + "\n")
	outputInfoFile.write("Unique Probes Condensed: " + str(uniqueNum) + "\n")
	for l in valueDict.keys():
		if len(valueDict[l]) > 1:
			outputInfoFile.write(l + "\t" + "\t".join(valueDict[l]) + "\n")
	outputInfoFile.close()

	for n in fasta.keys():
		outputFastaFile.write(n)
		outputFastaFile.write(fasta[n] + "\n")
	outputFastaFile.close()

def checkResults(finalFasta, checkBed, checkSortBed, originBed, uncovBed):
	finalResults = open(finalFasta, "r")
	checkBedFile = open(checkBed, "w")
	for line in finalResults:
		if line.startswith(">"):
			lines = line.split(">")[1].split("_")[0].split(":")
			chrom = lines[0]
			start_o = lines[1].split("-")[0]
			end = lines[3].split("-")
			start = int(start_o) + int(end[0])
			end = int(start_o) + int(end[1])
			checkBedFile.write(chrom + "\t" + str(start) + "\t" + str(end) + "\n")
	checkBedFile.close()
	cmd = """
		bedtools sort -i {checkBed} | bedtools merge -i - > {checkSortBed}
		bedtools sort -i {originBed} | bedtools merge -i - | \\
			bedtools subtract -a - -b {checkSortBed} > {uncovBed}
	""".format(checkBed=checkBed, checkSortBed=checkSortBed, originBed=originBed, uncovBed=uncovBed)
	os.system(cmd)

def main(inputFile, outputDir, length, thread, tiling):

	if not os.path.exists(outputDir):
		os.makedirs(outputDir)
	
	if "bed" in inputFile:
		sampleName = inputFile.split(".bed")[0]
		if "/" in sampleName:
			sampleName = sampleName.split("/")[1]
		outputFile = outputDir + "/" + sampleName + ".fasta"
		sort_merge(inputFile)
		fixLength(length)
		getFasta(outputFile)
	# elif "fasta" in inputFile:
	# 	sampleName = inputFile.split(".fasta")[0]
	# 	if "/" in sampleName:
	# 		sampleName = sampleName.split("/")[1]
	# 	outputFile = outputDir + "/" + sampleName + ".fasta"
	# 	shutil.copy(inputFile, outputFile)
	# elif "fa" in inputFile:
	# 	sampleName = inputFile.split(".fa")[0]
	# 	if "/" in sampleName:
	# 		sampleName = sampleName.split("/")[1]		
	# 	outputFile = outputDir + "/" + sampleName + ".fasta"
	# 	shutil.copy(inputFile, outputFile)
	else:
		print("请补充文件扩增名.bed")
		exit()

	outputFasta = outputDir + "/" + sampleName + "_" + str(tiling) + "X_Tiling.probes"
	UniqueFasta = outputFasta + ".final.fasta"
	outputInfo = outputDir + "/" + sampleName + "_" + str(tiling) + "X_Tiling.probes.report.txt"

	Design(outputFasta, outputFile, length, thread, tiling)
	os.remove(outputDir + "/" + sampleName + "_" + str(tiling) + "X_Tiling.probes.sqlite")
	rmDuplicates(outputFasta + ".fasta", UniqueFasta, outputInfo)

	checkBed = outputDir + "/" + sampleName + "_" + str(tiling) + "X_cov.bed"
	checkSortBed = outputDir + "/" + sampleName + "_" + str(tiling) + "X_cov.sort.bed"
	uncovBed = outputDir + "/" + sampleName + "_" + str(tiling) + "X_uncov.bed"
	checkResults(UniqueFasta, checkBed, checkSortBed, inputFile, uncovBed)

if __name__ == "__main__":
	parser = argparse.ArgumentParser(
		description="ProbeDesign  PZW",
		prog="ProbeDesign.py",
		usage="python3 ProbeDesign.py [-h] -i <bedfile> -f <fastaOutput> -o <outputFasta>",
		formatter_class=argparse.RawTextHelpFormatter
	)
	parser.add_argument("-v", "--version", action="version",
		version="Version 0.2 20200429")
	parser.add_argument("-i", "--input", type=str,
		help="可选择bed文件，必须含有.bed扩展名")
	parser.add_argument("-o", "--output", type=str,
		help="输出结果文件夹")
	parser.add_argument("-l", "--length", type=int,
		help="探针长度，默认为60bp", default=60)
	parser.add_argument("-t", "--thread", type=str,
		help="线程数，默认为8", default=8)
	parser.add_argument("-tile", "--tiling", type=int,
		help="叠瓦数，默认为3", default=3)
	if len(sys.argv[1:]) == 0:
		parser.print_help()
		parser.exit()
	args = parser.parse_args()
	main(inputFile=args.input, outputDir=args.output, length=args.length, thread=args.thread, tiling=args.tiling)