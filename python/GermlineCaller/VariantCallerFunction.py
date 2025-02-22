# pzw
# 20190620
# v0.4
# coding=utf-8

import os
import sys
import time
import argparse

def getAbsPath():
	now = os.path.abspath(os.path.dirname(sys.argv[0]))
	return now

def getTime():
	localtime = time.asctime(time.localtime(time.time()))
	return localtime

def sampleReady(sampleFile, rawdataPath):
	sampleID = open(sampleFile, "r")
	IDList = []
	for ID in sampleID:
		ID_fix = ID.split("\n")[0].split("\r")[0]
		IDList.append(ID_fix)

	fileList = os.listdir(rawdataPath)
	fileDict = {}
	for i in IDList:
		fileDict[i] = []
		for file in fileList:
			if (i + "_") in file or (i + ".") in file:
				if "R1" in file:
					fileDict[i].append(file)
				elif "R2" in file:
					fileDict[i].append(file)
				else:
					fileDict[i].append(file)
			else:
				continue

	return fileDict

def FastpfilterPair(input1, input2, output1, output2, jsonFile, htmlFile, thread):
	fastp = getAbsPath() + "/software/fastp"
	print("[filter start]", getTime())
	cmd = """
		{fastp} \\
			-i {input1} \\
			-I {input2} \\
			-o {output1} \\
			-O {output2} \\
			-j {jsonFile} \\
			-h {htmlFile} \\
			-w {thread} \\
			-q 5 -u 50 -n 15
	""".format(fastp=fastp, input1=input1, input2=input2, output1=output1, output2=output2, jsonFile=jsonFile, htmlFile=htmlFile, thread=thread)
	os.system(cmd)
	print("[filter done]", getTime())

def FastpfilterSingle(input1, output1, jsonFile, htmlFile, thread):
	fastp = getAbsPath() + "/software/fastp"
	print("[filter start]", getTime())
	cmd = """
		{fastp} \\
			-i {input1} \\
			-o {output1} \\
			-j {jsonFile} \\
			-h {htmlFile} \\
			-w {thread} \\
			-q 5 -u 50 -n 15
	""".format(fastp=fastp, input1=input1, output1=output1, jsonFile=jsonFile, htmlFile=htmlFile, thread=thread)
	os.system(cmd)
	print("[filter done]", getTime())


def mappingPair(input1, input2, output, thread):
	print("[mapping start]", getTime())
	bwa = getAbsPath() + "/software/bwa-0.7.17/bwa"
	samtools = getAbsPath() + "/software/samtools"
	reference = getAbsPath() + "/reference/hg19.fa"
	cmd = """
		{bwa} mem -t {thread} {reference} {input1} {input2} | {samtools} view -bSh - | \\
			{samtools} sort -@ {thread} - -o {output}
	""".format(bwa=bwa, thread=thread, reference=reference, input1=input1, input2=input2, output=output, samtools=samtools)
	os.system(cmd)
	print("[mapping done]", getTime())


def mappingSingle(input1, output, thread):
	print("[mapping start]", getTime())
	bwa = getAbsPath() + "/software/bwa-0.7.17/bwa"
	samtools = getAbsPath() + "/software/samtools"
	reference = getAbsPath() + "/reference/hg19.fa"
	cmd = """
		{bwa} mem -t {thread} {reference} {input1} | {samtools} view -bSh - | \\
			{samtools} sort -@ {thread} - -o {output}
	""".format(bwa=bwa, thread=thread, reference=reference, input1=input1, output=output, samtools=samtools)
	os.system(cmd)
	print("[mapping done]", getTime())


def makeBam(sample, outputDir):
	print("[Recalibrator start]", getTime())
	gatk = getAbsPath() + "/software/gatk-4.1.2.0/gatk"
	reference = getAbsPath() + "/reference/hg19.fa"
	indeldata1 = getAbsPath() + "/reference/Mills_and_1000G_gold_standard.indels.hg19.sites.vcf"
	indeldata2 = getAbsPath() + "/reference/1000G_phase1.indels.hg19.sites.vcf"
	cmd = """
		if [ ! -d "{outputDir}/bam" ]; then
			mkdir {outputDir}/bam
		fi
		{gatk} MarkDuplicates -I {outputDir}/temp/{sample}.bam -O {outputDir}/temp/{sample}.md.bam -M {outputDir}/temp/{sample}.md.txt
		{gatk} AddOrReplaceReadGroups -I {outputDir}/temp/{sample}.md.bam -O {outputDir}/temp/{sample}.ah.bam \\
			-LB lib1 -PL illumina -PU unit1 -SM {sample}
		{gatk} BuildBamIndex -I {outputDir}/temp/{sample}.ah.bam
		{gatk} BaseRecalibrator -I {outputDir}/temp/{sample}.ah.bam \\
			-R {reference} \\
			--known-sites {indeldata1} \\
			--known-sites {indeldata2} \\
			-O {outputDir}/temp/{sample}.recal.table
		{gatk} ApplyBQSR -R {reference} \\
			-I {outputDir}/temp/{sample}.ah.bam \\
			--bqsr {outputDir}/temp/{sample}.recal.table \\
			-O {outputDir}/bam/{sample}.final.bam
	""".format(gatk=gatk, reference=reference, indeldata1=indeldata1, indeldata2=indeldata2, sample=sample, outputDir=outputDir)
	os.system(cmd)
	print("[Recalibrator done]", getTime())


def convertBedToIntervalList(bed, intervalList):
	gatk = getAbsPath() + "/software/gatk-4.1.2.0/gatk"
	referDict = getAbsPath() + "/reference/hg19.dict"
	cmd = """
		{gatk} BedToIntervalList \\
			-SD {referDict} \\
			-I {bed} \\
			-O {intervalList}
	""".format(gatk=gatk, referDict=referDict, bed=bed, intervalList=intervalList)
	os.system(cmd)


def callsnp(inputBam, outputVcf, intervalList=""):
	print("[germline caller start]", getTime())
	gatk = getAbsPath() + "/software/gatk-4.1.2.0/gatk"
	reference = getAbsPath() + "/reference/hg19.fa"

	if intervalList == "":
		cmd = """
			{gatk} HaplotypeCaller \\
				-R {reference} \\
				-I {inputBam} \\
				-O {outputVcf}
		""".format(gatk=gatk, reference=reference, inputBam=inputBam, outputVcf=outputVcf)
	
	else:
		cmd = """
			{gatk} HaplotypeCaller \\
				-R {reference} \\
				-I {inputBam} \\
				-L {intervalList} \\
				-O {outputVcf}
		""".format(gatk=gatk, reference=reference, inputBam=inputBam, outputVcf=outputVcf, intervalList=intervalList)
	os.system(cmd)
	print("[germline caller done]", getTime())


def SomaticCall(inputBam, outputVcf, intervalList=""):
	print("[somatic caller start]", getTime())
	gatk = getAbsPath() + "/software/gatk-4.1.2.0/gatk"
	reference = getAbsPath() + "/reference/hg19.fa"
	gnomad = getAbsPath() + "/reference/af-only-gnomad.raw.sites.hg19.vcf.gz"
	
	if intervalList == "":
		cmd = """
			{gatk} Mutect2 -R {reference} \\
				-I {inputBam} \\
				--germline-resource {gnomad} \\
				-O {outputVcf}

		""".format(gatk=gatk, reference=reference, inputBam=inputBam, gnomad=gnomad, outputVcf=outputVcf)
	else:
		cmd = """
			{gatk} Mutect2 -R {reference} \\
				-I {inputBam} \\
				--germline-resource {gnomad} \\
				-L {intervalList} \\
				-O {outputVcf}

		""".format(gatk=gatk, reference=reference, inputBam=inputBam, gnomad=gnomad, outputVcf=outputVcf, intervalList=intervalList)	

	os.system(cmd)
	print("[somatic caller done", getTime())


def annotate(outputDir, sample, thread):
	print("[annotate start]", getTime())
	annovarDir = getAbsPath() + "/software/annovar2018"
	humandb = annovarDir + "/humandb/"
	cmd = """
		{annovarDir}/convert2annovar.pl -format vcf4 {outputDir}/vcf/{sample}.vcf > {outputDir}/annotation/avinput
		{annovarDir}/table_annovar.pl {outputDir}/annotation/avinput {humandb} -buildver hg19 \\
			-out {outputDir}/annotation/{sample} -remove \\
			-protocol refGene,cytoBand,avsnp150,1000g2015aug_all,1000g2015aug_eas,exac03,clinvar_20190305,cosmic70,dbnsfp35a \\
			-operation g,r,f,f,f,f,f,f,f \\
			-nastring . -thread {thread} -otherinfo
	""".format(outputDir=outputDir, sample=sample, annovarDir=annovarDir, thread=thread, humandb=humandb)
	os.system(cmd)
	print("[annotate done]", getTime())


def main(inputDir, outputDir, sampleList, bed, threads, method):
	process = sampleReady(sampleList, inputDir)
	for sample in process:
		print(sample, "begin: ", getTime())
		if len(process[sample]) == 2:
			read1 = inputDir + "/" + process[sample][0]
			read2 = inputDir + "/" + process[sample][1]
			cmd = """
				if [ ! -d "{outputDir}" ]; then
					mkdir {outputDir}
				fi

				if [ ! -d "{outputDir}/QC" ]; then
					mkdir {outputDir}/QC
				fi

				if [ ! -d "{outputDir}/cleandata" ]; then
					mkdir {outputDir}/cleandata
				fi

			""".format(outputDir=outputDir)
			os.system(cmd)


			FastpfilterPair(
				read1, read2,
				outputDir + "/cleandata/{sample}.R1.fq.gz".format(sample=sample),
				outputDir + "/cleandata/{sample}.R2.fq.gz".format(sample=sample),
				outputDir + "/QC/{sample}.json".format(sample=sample),
				outputDir + "/QC/{sample}.html".format(sample=sample),
				threads
			)


			cmd = """
				if [ ! -d "{outputDir}/temp" ]; then
					mkdir {outputDir}/temp
				fi
			""".format(outputDir=outputDir)
			os.system(cmd)

			mappingPair(
				outputDir + "/cleandata/{sample}.R1.fq.gz".format(sample=sample),
				outputDir + "/cleandata/{sample}.R2.fq.gz".format(sample=sample),
				outputDir + "/temp/{sample}.bam".format(sample=sample),
				threads
			)

		elif len(process[sample]) == 1:
			print("only find 1 rawdata, guess its single end.")
			read = inputDir + "/" + process[sample][0]
			cmd = """
				if [ ! -d "{outputDir}" ]; then
					mkdir {outputDir}
				fi

				if [ ! -d "{outputDir}/QC" ]; then
					mkdir {outputDir}/QC
				fi

				if [ ! -d "{outputDir}/cleandata" ]; then
					mkdir {outputDir}/cleandata
				fi

			""".format(outputDir=outputDir)
			os.system(cmd)

			FastpfilterSingle(
				read,
				outputDir + "/cleandata/{sample}.fq.gz".format(sample=sample),
				outputDir + "/QC/{sample}.json".format(sample=sample),
				outputDir + "/QC/{sample}.html".format(sample=sample),
				threads
			)


			cmd = """
				if [ ! -d "{outputDir}/temp" ]; then
					mkdir {outputDir}/temp
				fi
			""".format(outputDir=outputDir)
			os.system(cmd)

			mappingSingle(
				outputDir + "/cleandata/{sample}.fq.gz".format(sample=sample),
				outputDir + "/temp/{sample}.bam".format(sample=sample),
				threads
			)
		elif len(process[sample]) == 0:
			print("can not find sample, please check sample ID or rawdata")
			continue

		else:
			print("find more than 2 rawdata, skip this sample")
			continue

		makeBam(sample, outputDir)

		cmd = """
			if [ ! -d "{outputDir}/vcf" ]; then
				mkdir {outputDir}/vcf
			fi
		""".format(outputDir=outputDir)
		os.system(cmd)

		if bed != "":
			convertBedToIntervalList(bed, "{outputDir}/temp/temp.list".format(outputDir=outputDir))
			if method == "germline":
				callsnp(
					outputDir + "/bam/{sample}.final.bam".format(sample=sample),
					outputDir + "/vcf/{sample}.vcf".format(sample=sample),
					outputDir + "/temp/temp.list"
				)
			elif method == "somatic":
				SomaticCall(
					outputDir + "/bam/{sample}.final.bam".format(sample=sample),
					outputDir + "/vcf/{sample}.vcf".format(sample=sample),
					outputDir + "/temp/temp.list"
				)
			else:
				print("ERROR, please select a method")
		else:
			if method == "germline":
				callsnp(
					outputDir + "/bam/{sample}.final.bam".format(sample=sample),
					outputDir + "/vcf/{sample}.vcf".format(sample=sample)
				)
			elif method == "somatic":
				SomaticCall(
					outputDir + "/bam/{sample}.final.bam".format(sample=sample),
					outputDir + "/vcf/{sample}.vcf".format(sample=sample)
				)
			else:
				print("ERROR, please select a method")				

		cmd = """
			if [ ! -d "{outputDir}/annotation" ]; then
				mkdir {outputDir}/annotation
			fi
			rm -rf {outputDir}/temp
		""".format(outputDir=outputDir)
		os.system(cmd)

		annotate(outputDir, sample, threads)
		cmd = """
			rm {outputDir}/annotation/avinput
		""".format(outputDir=outputDir)
		os.system(cmd)
		print(sample, "finish: ", getTime())

if __name__ == "__main__":
	parser = argparse.ArgumentParser(description="Variant Caller Pipeline",
		prog="VariantCaller",
		usage="VariantCaller -i <input dir> -o <output dir> -s <sample list> -t <threads> -m <method>")
	group = parser.add_mutually_exclusive_group()
	parser.add_argument("-v", "--version", action="version",
		version="Version 1.0 20190610")
	parser.add_argument("-i", "--input", type=str,
		help="Input the directory of the fq.gz file")
	parser.add_argument("-o", "--output", type=str,
		help="the output directory")
	parser.add_argument("-s", "--samplelist", type=str,
		help="sample list, one id in each row")
	parser.add_argument("-m", "--method", type=str, default="",
		help="select variant type, 'somatic' or 'germline'")
	parser.add_argument("-b", "--bed", type=str, default="",
		help="input bed file, optional")
	parser.add_argument("-t", "--threads", type=int, default=4,
		help="threads, default=4")
	if len(sys.argv[1:]) == 0:
		parser.print_help()
		parser.exit()
	args = parser.parse_args()
	main(inputDir=args.input, outputDir=args.output, sampleList=args.samplelist, bed=args.bed, threads=args.threads, method=args.method)
