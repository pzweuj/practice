# coding=utf-8
# pzw
# 20191127

import sys
import os
import function.WordWriter3 as ww
import time
import argparse
from collections import Counter

def getAbsPath():
	now = os.path.abspath(os.path.dirname(sys.argv[0]))
	return now

def sampleinfo(samplefile):
	sample = open(samplefile, "r")
	sampleDict = {}
	for line in sample:
		if line.startswith("受检人姓名"):
			continue
		else:
			sampleString = line.split("\n")[0].split(",")
			sampleDict[sampleString[3]] = {
				"name": sampleString[0],
				"gender": sampleString[1],
				"age": sampleString[2],
				"id": sampleString[3],
				"type": sampleString[4],
				"sampleTime": sampleString[5],
				"hospital_id": sampleString[6],
				"company": sampleString[7],
				"hospital_part": sampleString[8],
				"doctor": sampleString[9],
				"phone": sampleString[10],
				"diagnosis": sampleString[11],
				"history": sampleString[12],
				"drug": sampleString[13]
			}
	sample.close()
	return sampleDict


def filter(fastq1, fastq2, ids, outputDir):
	cmd = """
		fastp -i {fastq1} -I {fastq2} \\
			-o {outputDir}/{ids}_1.fastq.gz \\
			-O {outputDir}/{ids}_2.fastq.gz \\
			-j {outputDir}/{ids}.json \\
			-h {outputDir}/{ids}.html \\
			-w 8

	""".format(fastq1=fastq1, fastq2=fastq2, ids=ids, outputDir=outputDir)
	os.system(cmd)

def filter2(fastq1, ids, outputDir):
	cmd = """
		fastp -i {fastq1} \\
			-o {outputDir}/{ids}_1.fastq.gz \\
			-j {outputDir}/{ids}.json \\
			-h {outputDir}/{ids}.html \\
			-w 8
	""".format(fastq1=fastq1, ids=ids, outputDir=outputDir)
	os.system(cmd)


def bowtie2mapping(fastq1, fastq2, ids, outputDir):
	now = getAbsPath()
	reference = now + "/reference/HBVRT"
	cmd = """
		bowtie2 -x {reference} -1 {fastq1} -2 {fastq2} \\
			-S {outputDir}/{ids}.sam --very-sensitive --end-to-end -p 8
		samtools view {outputDir}/{ids}.sam -F 12 -bSh | samtools sort -@ 8 - -o {outputDir}/{ids}.bam
		rm {outputDir}/{ids}.sam
	""".format(fastq1=fastq1, fastq2=fastq2, ids=ids, outputDir=outputDir, reference=reference)
	os.system(cmd)

def bowtie2mapping2(fastq1, ids, outputDir):
	now = getAbsPath()
	reference = now + "/reference/HBVRT"
	cmd = """
		bowtie2 -x {reference} -U {fastq1} \\
			-S {outputDir}/{ids}.sam --very-sensitive --end-to-end -p 8

		samtools view {outputDir}/{ids}.sam -F 12 -bSh | samtools sort -@ 8 - -o {outputDir}/{ids}.bam
		rm {outputDir}/{ids}.sam
	""".format(fastq1=fastq1, ids=ids, outputDir=outputDir, reference=reference)
	os.system(cmd)

def blastResults(inputBam, outputDir, ids, cutoff):
	now = getAbsPath()
	blastdb = now + "/reference/HBV48.fasta"
	typing = now + "/function/HBV_typing2.py"

	# reads阈值，只有超过这个reads数才认为是检出阳性
	# cutoff = "20"

	cmd = """
		bedtools bamtofastq -i {inputBam} -fq {outputDir}/{ids}.fastq
		cat {outputDir}/{ids}.fastq | paste - - - - | sed 's/^@/>/g' | cut -f1-2 | tr '\\t' '\\n' > {outputDir}/{ids}.fasta
		rm {outputDir}/{ids}.fastq
	""".format(inputBam=inputBam, outputDir=outputDir, ids=ids)
	os.system(cmd)

	queryTemp = []
	fastaFix = open(outputDir + "/" + ids + ".fasta", "r")
	for line in fastaFix:
		if line.startswith(">"):
			continue
		else:
			queryTemp.append(line)
	queryDict = Counter(queryTemp)
	fastaFixed = open(outputDir + "/" + ids + ".fix.fasta", "w")
	uniqueDict = {}
	n = 1
	for i in queryDict.keys():
		counts = queryDict[i]
		readID = ">uniqueID_" + str(n) + "_counts_" + str(counts)
		uniqueDict[readID] = counts
		n += 1
		fastaFixed.write(readID + "\n")
		fastaFixed.write(i)
	fastaFixed.close()
	fastaFix.close()

	cmd2 = """
		blastn -query {outputDir}/{ids}.fix.fasta -db {blastdb} -out {outputDir}/{ids}.blast.txt -num_threads 8
		rm {outputDir}/{ids}.fasta {outputDir}/{ids}.fix.fasta
		python2 {typing} -i {outputDir}/{ids}.blast.txt -c {cutoff} > {outputDir}/{ids}.typingResults.txt
		rm {outputDir}/{ids}.blast.txt
	""".format(inputBam=inputBam, outputDir=outputDir, ids=ids, blastdb=blastdb, typing=typing, cutoff=cutoff)
	os.system(cmd2)

def drugVcf(inputBam, outputDir, ids, cutoff, skipIns):
	now = getAbsPath()
	reference = now + "/reference/HBVRT.fasta"
	drug = now + "/function/HBV_drug2.py"
	drugDB = now + "/reference/HBV_drug2.txt"

	# 深度阈值，只有超过这个深度才认为是真实突变
	# cutoff = "200"
	if skipIns:
		cmd = """
			bcftools mpileup -f {reference} {inputBam} | bcftools call -mv -O v -o {outputDir}/{ids}.vcf
			python2 {drug} -i {outputDir}/{ids}.vcf -c {cutoff} -r {reference} -d {drugDB} -s True > {outputDir}/{ids}.drug.txt
			python2 {drug} -i {outputDir}/{ids}.vcf -c {cutoff} -r {reference} -d {drugDB} -m list -s True > {outputDir}/{ids}.drug.list.txt
			python2 {drug} -i {outputDir}/{ids}.vcf -c {cutoff} -r {reference} -d {drugDB} -m seq -s True > {outputDir}/{ids}.S.txt
		""".format(reference=reference, inputBam=inputBam, outputDir=outputDir, ids=ids, drug=drug, drugDB=drugDB, cutoff=cutoff)
		os.system(cmd)
	else:
		cmd = """
			bcftools mpileup -f {reference} {inputBam} | bcftools call -mv -O v -o {outputDir}/{ids}.vcf
			python2 {drug} -i {outputDir}/{ids}.vcf -c {cutoff} -r {reference} -d {drugDB} > {outputDir}/{ids}.drug.txt
			python2 {drug} -i {outputDir}/{ids}.vcf -c {cutoff} -r {reference} -d {drugDB} -m list > {outputDir}/{ids}.drug.list.txt
			python2 {drug} -i {outputDir}/{ids}.vcf -c {cutoff} -r {reference} -d {drugDB} -m seq > {outputDir}/{ids}.S.txt
		""".format(reference=reference, inputBam=inputBam, outputDir=outputDir, ids=ids, drug=drug, drugDB=drugDB, cutoff=cutoff)
		os.system(cmd)


def fillReportDict(samplefile, ids, drugFile, drugListFile, typeFile):
	sampleDicts = sampleinfo(samplefile)
	fillDict = {}
	fillDict["#[TBS-NAME]#"] = sampleDicts[ids]["name"]
	fillDict["#[TBS-AGE]#"] = sampleDicts[ids]["age"]
	fillDict["#[TBS-SAMPLEID]#"] = ids
	fillDict["#[TBS-FROM]#"] = sampleDicts[ids]["company"]
	fillDict["#[TBS-DOCTOR]#"] = sampleDicts[ids]["doctor"]
	fillDict["#[TBS-DIAGNOISE]#"] = sampleDicts[ids]["diagnosis"]
	fillDict["#[TBS-ILLHISTORY]#"] = sampleDicts[ids]["history"]
	fillDict["#[TBS-DRUGHISTORY]#"] = sampleDicts[ids]["drug"]
	fillDict["#[TBS-GENDER]#"] = sampleDicts[ids]["gender"]
	fillDict["#[TBS-HOSPITALID]#"] = sampleDicts[ids]["hospital_id"]
	fillDict["#[TBS-COLLECTDATE]#"] = sampleDicts[ids]["sampleTime"]
	fillDict["#[TBS-HOSPITALPART]#"] = sampleDicts[ids]["hospital_part"]
	fillDict["#[TBS-PHONE]#"] = sampleDicts[ids]["phone"]

	typingFile = open(typeFile, "r")
	typesList = []
	for line in typingFile:
		types = line.split("\n")[0].split("\t")[0]
		if "未" in types:
			continue
		else:
			typesList.append(types)
	fillDict["#[TYPINGAMOUNT]#"] = str(len(typesList))
	if len(typesList) == 0:
		fillDict["#[TYPE]#"] = "本次检测未检测到明确分型"
	else:
		fillDict["#[TYPE]#"] = "，".join(typesList)
	typingFile.close()

	fillDict["#[TABLE-typing]#"] = typeFile

	fillDict["#[HEADER-name]#"] = sampleDicts[ids]["name"]
	localtime = time.localtime()
	fillDict["#[HEADER-date]#"] = time.strftime("%Y年%m月%d日", localtime)

	drug = open(drugFile, "r")
	for line in drug:
		if line == "\n":
			continue
		elif "未发现" in line:
			fillDict["#[BOOLDRUG]#"] = "未"
			fillDict["#[DRUGAMOUNT]#"] = ""
			fillDict["#[FILLSUPP]#"] = ""
			fillDict["#[DRUGS]#"] = ""
		else:
			fillDict["#[BOOLDRUG]#"] = ""
			drugs = line.split("\n")[0]
			fillDict["#[FILLSUPP]#"] = "，可能对以下药物耐药："
			fillDict["#[DRUGS]#"] = drugs

			drugAAchangeAmount = []
			drugAAchangeAmountFile = open(drugListFile, "r")
			for l in drugAAchangeAmountFile:
				if "未突变" in l:
					continue
				else:
					drugAAchangeAmount.append(l)
			drugAAchangeAmount = list(set(drugAAchangeAmount))
			fillDict["#[DRUGAMOUNT]#"] = str(len(drugAAchangeAmount)) + "个"
	drug.close()

	fillDict["#[TABLE-drug]#"] = drugListFile
	return fillDict

def main(sampleinfoFile, ids, fastq1, fastq2, outputDir, vcfCutoff, readCutOff, skipIns):


	sampleDictionary = sampleinfo(sampleinfoFile)
	print (ids + " [start]")
	cmd = """
		mkdir -p {outputDir}/{ids}
		mkdir {outputDir}/{ids}/report
		mkdir {outputDir}/{ids}/fillUp
		mkdir {outputDir}/{ids}/cleandata
		mkdir {outputDir}/{ids}/bam
	""".format(outputDir=outputDir, ids=ids)
	os.system(cmd)

	print (ids + " [filter]")
	if fastq2 == "single":
		filter2(fastq1, ids, outputDir + "/" + ids + "/cleandata")
	else:
		filter(fastq1, fastq2, ids, outputDir + "/" + ids + "/cleandata")

	print (ids + " [mapping]")
	cleanfq1 = outputDir + "/" + ids + "/cleandata/" + ids + "_1.fastq.gz"
	cleanfq2 = outputDir + "/" + ids + "/cleandata/" + ids + "_2.fastq.gz"
	if fastq2 == "single":
		bowtie2mapping2(cleanfq1, ids, outputDir + "/" + ids + "/bam")
	else:
		bowtie2mapping(cleanfq1, cleanfq2, ids, outputDir + "/" + ids + "/bam")

	print (ids + " [blast]")
	fillTemp = outputDir + "/" + ids + "/fillUp"
	bamfile = outputDir + "/" + ids + "/bam/" + ids + ".bam"
	blastResults(bamfile, fillTemp, ids, readCutOff)

	print (ids + " [calling]")
	drugVcf(bamfile, fillTemp, ids, vcfCutoff, skipIns)

	print (ids + " [analysis]")
	drugFile = fillTemp + "/" + ids + ".drug.txt"
	drugListFile = fillTemp + "/" + ids + ".drug.list.txt"
	typeFile = fillTemp + "/" + ids + ".typingResults.txt"
	reportDict = fillReportDict(sampleinfoFile, ids, drugFile, drugListFile, typeFile)

	print (ids + " [reporting]")
	doxtemp = getAbsPath() + "/template/HBV_20191114.docx"
	outputdocx = outputDir + "/" + ids + "/report/" + ids + "-HBV-" + reportDict["#[TBS-NAME]#"] + "-" + reportDict["#[HEADER-date]#"] + ".docx"
	ww.WordWriter(doxtemp, outputdocx, reportDict)
	print (ids + " [finish]")
		
if __name__ == "__main__":
	parser = argparse.ArgumentParser(
		description="HBV typing and drug analysis  PZW@Genephar",
		prog="HBVpro.py",
		usage="python3 HBVpro.py [-h] -s <sampleinfo> -i <sampleID> -f1 <fastq1> -f2 <fastq2> -o <outputDir>",
		formatter_class=argparse.RawTextHelpFormatter
	)
	parser.add_argument("-v", "--version", action="version",
		version="Version 2.0 20191127")
	parser.add_argument("-s", "--sampleinfo", type=str,
		help="Input the sampleinfo file", default=getAbsPath() + "/reference/HBV.txt")
	parser.add_argument("-i", "--ids", type=str,
		help="Input the sample ID", default="H00001")
	parser.add_argument("-f1", "--fastq1", type=str,
		help="Input _1_fastq.gz file")
	parser.add_argument("-f2", "--fastq2", type=str,
		help="Input _2_fastq.gz file", default="single")
	parser.add_argument("-o", "--output", type=str,
		help="specify a output directory")
	parser.add_argument("-vcf", "--vcfCutoff", type=str,
		help="vcf cutoff, default=20", default="20")
	parser.add_argument("-r", "--readsCutoff", type=str,
		help="reads cutoff, default=20", default="20")
	parser.add_argument("-skip", "--skipins", type=bool,
		help="skip indel, default=False", default=False)
	if len(sys.argv[1:]) == 0:
		parser.print_help()
		parser.exit()
	args = parser.parse_args()
	main(sampleinfoFile=args.sampleinfo, ids=args.ids, fastq1=args.fastq1, fastq2=args.fastq2, outputDir=args.output, vcfCutoff=args.vcfCutoff, readCutOff=args.readsCutoff, skipIns=args.skipins)