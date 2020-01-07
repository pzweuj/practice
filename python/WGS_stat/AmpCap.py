# coding:utf-8
# pzw
# 20200107
# AmpCap
# amplicon capture analysis

import os
import sys
import json
import subprocess
import argparse

def fastp(inputDir, outputDir, sample, thread, adapter="auto", trim="0"):
	folder = os.path.exists(outputDir + "/QC")
	if not folder:
		os.makedirs(outputDir + "/QC")

	folder = os.path.exists(outputDir + "/cleandata")
	if not folder:
		os.makedirs(outputDir + "/cleandata")

	fastp = "fastp"
	cmd = """
		{fastp} -i {inputDir}/{sample}_R1_001.fastq.gz \\
			-I {inputDir}/{sample}_R2_001.fastq.gz \\
			-o {outputDir}/cleandata/{sample}_1.fq.gz \\
			-O {outputDir}/cleandata/{sample}_2.fq.gz \\
			-w {thread} -h {outputDir}/QC/{sample}.html -j {outputDir}/QC/{sample}.json \\
			-a {adapter} -f {trim}
	""".format(fastp=fastp, inputDir=inputDir, sample=sample, outputDir=outputDir, thread=thread, adapter=adapter, trim=trim)
	
	os.system(cmd)


def bwa(outputDir, sample, thread):
	folder = os.path.exists(outputDir + "/bam")
	if not folder:
		os.makedirs(outputDir + "/bam")
	
	bwa = "bwa"
	samtools = "samtools"
	hg19 = "hg19.fa"
	cmd = """
		{bwa} mem -t {thread} {hg19} \\
			{outputDir}/cleandata/{sample}_1.fq.gz \\
			{outputDir}/cleandata/{sample}_2.fq.gz \\
			| {samtools} view -bSh - | {samtools} sort -@ {thread} - -o {outputDir}/bam/{sample}.bam
	""".format(bwa=bwa, thread=thread, hg19=hg19, outputDir=outputDir, sample=sample, samtools=samtools)

	os.system(cmd)

def bedtools(outputDir, sample):
	folder = os.path.exists(outputDir + "/covStat")
	if not folder:
		os.makedirs(outputDir + "/covStat")

	bedtools = "bedtools"
	cmd = """
		{bedtools} genomecov -ibam {outputDir}/bam/{sample}.bam -bga > {outputDir}/covStat/{sample}.cov.txt
	""".format(bedtools=bedtools, outputDir=outputDir, sample=sample)

	os.system(cmd)


def LocationCov(outputDir, sample):
	covFile = open(outputDir + "/covStat/" + sample + ".cov.txt", "r")
	covDict = {}
	for line in covFile:
		lines = line.split("\n")[0].split("\t")
		if int(lines[3]) == 0:
			continue
		else:
			for i in range(int(lines[1])+1, int(lines[2])+1):
				covDict[lines[0] + ":" + str(i)] = lines[3]
	covFile.close()
	return covDict

def getMappingReads(outputDir, sample):
	samtools = "samtools"
	cmd = """
		{samtools} view {outputDir}/bam/{sample}.bam -F 2052 | wc -l
	""".format(samtools=samtools, outputDir=outputDir, sample=sample)
	results = subprocess.check_output(cmd, shell=True)
	return str(results).split("b'")[1].split("\\n")[0]


def bamAnalysis(outputDir, sample, bedFile):

	coverageDict = LocationCov(outputDir, sample)
	bed = open(bedFile, "r")

	outputFile = open(outputDir + "/covStat/" + sample + ".results.txt", "w")
	jsonFile = json.load(open(outputDir + "/QC/" + sample + ".json", "r"))
	rawReads = jsonFile["summary"]["before_filtering"]["total_reads"]
	rawBases = jsonFile["summary"]["before_filtering"]["total_bases"]
	cleanReads = jsonFile["summary"]["after_filtering"]["total_reads"]
	cleanBases = jsonFile["summary"]["after_filtering"]["total_bases"]
	mappingReads = getMappingReads(outputDir, sample)

	outputFile.write("# " + sample + "\n")
	outputFile.write("# Panel\t" + bedFile + "\n")
	outputFile.write("rawReads\t" + str(rawReads) + "\n")
	outputFile.write("rawBases\t" + str(rawBases) + "\n")
	outputFile.write("cleanReads\t" + str(cleanReads) + "\n")
	outputFile.write("cleanBases\t" + str(cleanBases) + "\n")
	outputFile.write("mappingReads\t" + mappingReads + "\n\n\n")


	outputFile.write("# " + "\t".join(["AmpID", "chrom", "start", "end", "length", "averageDepth", "coverage", "是否达到均一标准", "比例\n"]))

	depthSum = 0
	basesSum = 0
	lengthSum = 0

	outputLine = []

	for line in bed:
		lines = line.split("\t")
		chrom = lines[0]
		start = lines[1]
		end = lines[2]
		name = lines[3]

		covBases = 0
		baseNotCover = 0

		for i in range(int(start), int(end)+1):
			if not (chrom + ":" + str(i)) in coverageDict.keys():
				s = 0
				baseNotCover += 1
			else:
				s = int(coverageDict[chrom + ":" + str(i)])
			covBases += s

		length = int(end) - int(start) + 1
		averageDepth = "%.2f" % (float(covBases) / float(length))
		depthSum += float(covBases) / float(length)
		basesSum += covBases
		lengthSum += length

		covPercent = "%.2f" % (float((length - baseNotCover)) / float(length) * 100) + "%"

		outputLine.append([name.split("\n")[0], chrom, start, end, str(length), averageDepth, covPercent])

	for i in outputLine:
		if float(i[5]) >= (float(basesSum) / float(lengthSum) * 0.2):
			uniform = "yes"
		else:
			uniform = "no"
		if depthSum == 0:
			ReadPercent = "0"
		else:
			ReadPercent = str(float(i[5]) / depthSum)
		outputFile.write("\t".join(i) + "\t" + uniform + "\t" + ReadPercent + "\n")

	outputFile.close()
	bed.close()

def main(inputDir, outputDir, sample, thread, adapter, trim, bedFile):
	print("[质控中]")
	fastp(inputDir, outputDir, sample, thread, adapter, trim)
	print("[比对中]")
	bwa(outputDir, sample, thread)
	print("[提取信息中]")
	bedtools(outputDir, sample)
	print("[获得结果中]")
	bamAnalysis(outputDir, sample, bedFile)
	print("[完成]")


if __name__ == "__main__":
	parser = argparse.ArgumentParser(
		description="Amplicon Capture Analysis  PZW@Genephar",
		prog="AmpCap.py",
		usage="python3 AmpCap.py [-h] -s <sampleID> -i <rawdataDir> -o <outputDir> -b <bedFile>",
		formatter_class=argparse.RawTextHelpFormatter
	)
	parser.add_argument("-v", "--version", action="version",
		version="Version 0.2 20200107")
	parser.add_argument("-i", "--rawdata", type=str,
		help="Input the rawdata path")
	parser.add_argument("-o", "--output", type=str,
		help="specify output directory")
	parser.add_argument("-s", "--sample", type=str,
		help="specify a sample ID.")
	parser.add_argument("-b", "--bedfile", type=str,
		help="specify a bed file. format: chrom\tstart\tend\tname")
	parser.add_argument("-t", "--thread", type=str,
		help="threads to use, option, default=8", default="8")
	parser.add_argument("-a", "--adapter", type=str,
		help="adapter seq, option, default=auto detect", default="auto")
	parser.add_argument("-trim", "--trim", type=str,
		help="trim reads bases, option, default=0", default="0")	

	if len(sys.argv[1:]) == 0:
		parser.print_help()
		parser.exit()
	args = parser.parse_args()
	main(inputDir=args.rawdata, outputDir=args.output, sample=args.sample, thread=args.thread, adapter=args.adapter, trim=args.trim, bedFile=args.bedfile)