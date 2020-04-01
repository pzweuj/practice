# coding=utf-8
# pzw
# 20200330

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
	fastp = "fastp"
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

def mappingPair(input1, input2, output, thread):
	print("[mapping start]", getTime())
	bwa = "bwa"
	samtools = "samtools"
	reference = "hg19.fa"
	cmd = """
		{bwa} mem -t {thread} {reference} {input1} {input2} | {samtools} view -bSh - | \\
			{samtools} sort -@ {thread} - -o {output}
	""".format(bwa=bwa, thread=thread, reference=reference, input1=input1, input2=input2, output=output, samtools=samtools)
	os.system(cmd)
	print("[mapping done]", getTime())

def makeBam(sample, outputDir):
	print("[Recalibrator start]", getTime())
	gatk = "gatk"
	reference = "hg19.fa"
	indeldata1 = "Mills_and_1000G_gold_standard.indels.hg19.sites.vcf"
	indeldata2 = "1000G_phase1.indels.hg19.sites.vcf"
	cmd = """
		if [ ! -d "{outputDir}/bam" ]; then
			mkdir {outputDir}/bam
		fi
		{gatk} AddOrReplaceReadGroups -I {outputDir}/temp/{sample}.bam -O {outputDir}/temp/{sample}.ah.bam \\
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

def SomaticCall(inputBam, outputVcf, outputFixVcf, threads):
	print("[somatic caller start]", getTime())
	gatk = "gatk"
	reference = "hg19.fa"
	gnomad = "af-only-gnomad.raw.sites.hg19.vcf.gz"
	
	cmd = """
		{gatk} Mutect2 -R {reference} \\
			-I {inputBam} \\
			--germline-resource {gnomad} \\
			--native-pair-hmm-threads {threads} \\
			-O {outputVcf}
		{gatk} AnnotateVcfWithBamDepth -V {outputVcf} \\
			-I {inputBam} -O {outputFixVcf}

	""".format(threads=threads, gatk=gatk, reference=reference, inputBam=inputBam, gnomad=gnomad, outputVcf=outputVcf, outputFixVcf=outputFixVcf)

	os.system(cmd)
	print("[somatic caller done", getTime())

def annotate(outputDir, sample, thread):
	print("[annotate start]", getTime())
	annovarDir = "annovar2019"
	humandb = annovarDir + "/humandb/"
	cmd = """
		{annovarDir}/convert2annovar.pl -format vcf4 {outputDir}/vcf/{sample}.vcf -includeinfo > {outputDir}/annotation/avinput
		{annovarDir}/table_annovar.pl {outputDir}/annotation/avinput {humandb} -buildver hg19 \\
			-out {outputDir}/annotation/{sample} -remove \\
			-protocol refGene,cytoBand,avsnp150 \\
			-operation g,r,f \\
			-nastring . -thread {thread} -otherinfo
	""".format(outputDir=outputDir, sample=sample, annovarDir=annovarDir, thread=thread, humandb=humandb)
	os.system(cmd)
	print("[annotate done]", getTime())

def main(inputDir, outputDir, sampleList, threads):
	process = sampleReady(sampleList, inputDir)
	for sample in process:
		print(sample, "begin: ", getTime())
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
		makeBam(sample, outputDir)

		cmd = """
			if [ ! -d "{outputDir}/vcf" ]; then
				mkdir {outputDir}/vcf
			fi
		""".format(outputDir=outputDir)
		os.system(cmd)


		SomaticCall(
			outputDir + "/bam/{sample}.final.bam".format(sample=sample),
			outputDir + "/temp/{sample}.vcf".format(sample=sample),
			outputDir + "/vcf/{sample}.vcf".format(sample=sample),
			threads
		)		

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
	parser = argparse.ArgumentParser(description="Variant Caller Pipeline Thyroid Special Version",
		prog="thyroid_variants_call.py",
		usage="python3 thyroid_variants_call.py -i <input dir> -o <output dir> -s <sample list> -t <threads>")
	group = parser.add_mutually_exclusive_group()
	parser.add_argument("-v", "--version", action="version",
		version="Version 0.1 20200330")
	parser.add_argument("-i", "--input", type=str,
		help="Input the directory of the fq.gz file")
	parser.add_argument("-o", "--output", type=str,
		help="the output directory")
	parser.add_argument("-s", "--samplelist", type=str,
		help="sample list, one id in each row")
	parser.add_argument("-t", "--threads", type=int, default=4,
		help="threads, default=4")
	if len(sys.argv[1:]) == 0:
		parser.print_help()
		parser.exit()
	args = parser.parse_args()
	main(inputDir=args.input, outputDir=args.output, sampleList=args.samplelist, threads=args.threads)
