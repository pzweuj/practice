# coding=utf-8
# pzw
# 20200317

import os
import sys
import argparse
import commands

# 原始数据质控
def filter(sample, rawdataDir, outputDir):
	fastp = "fastp"
	cmd = """
		if [ ! -d "{outputDir}" ]; then
			mkdir {outputDir}
		fi
		mkdir {outputDir}/cleandata
		mkdir {outputDir}/QC
		{fastp} -i {rawdataDir}/{sample}_R1_001.fastq.gz \\
			-I {rawdataDir}/{sample}_R2_001.fastq.gz \\
			-o {outputDir}/cleandata/{sample}_1.fq.gz \\
			-O {outputDir}/cleandata/{sample}_2.fq.gz \\
			-w 8 -j {outputDir}/QC/{sample}.json \\
			-h {outputDir}/QC/{sample}.html
	""".format(fastp=fastp, sample=sample, rawdataDir=rawdataDir, outputDir=outputDir)
	os.system(cmd)

# 比对
def mapping(sample, rawdataDir, outputDir, QC):
	reference = "hg19.fa"
	bwa = "bwa"
	if QC == "yes":
		cmd = """
			mkdir {outputDir}/bam
			{bwa} mem -t 8 {reference} {outputDir}/cleandata/{sample}_1.fq.gz \\
				{outputDir}/cleandata/{sample}_2.fq.gz \\
				| samtools view -bSh - | samtools sort -@ 8 - -o {outputDir}/bam/{sample}.bam
			samtools index {outputDir}/bam/{sample}.bam
		""".format(bwa=bwa, reference=reference, outputDir=outputDir, sample=sample)
	else:
		cmd = """
			mkdir {outputDir}/bam
			{bwa} mem -t 8 {reference} {rawdataDir}/{sample}_R1_001.fastq.gz \\
				{rawdataDir}/{sample}_R2_001.fastq.gz \\
				| samtools view -bSh - | samtools sort -@ 8 - -o {outputDir}/bam/{sample}.bam
			samtools index {outputDir}/bam/{sample}.bam
		""".format(bwa=bwa, reference=reference, rawdataDir=rawdataDir, outputDir=outputDir, sample=sample)
	os.system(cmd)

# 捕获分析
def amplicon(sample, outputDir, panel):
	bamdst = "bamdst"
	bedFile = "hg19.msi.bed"
	if panel == "2":
		bedFile = "hg19.msisensor.test.bed"
	bedtools = "bedtools"
	cmd = """
		mkdir {outputDir}/amplicon
		mkdir {outputDir}/amplicon/temp
		{bamdst} -p {bedFile} -o {outputDir}/amplicon/temp {outputDir}/bam/{sample}.bam
		zcat {outputDir}/amplicon/temp/region.tsv.gz \\
			| {bedtools} intersect -a - -b {bedFile} -wb -header \\
			| awk 'BEGIN{{OFS="\t"}} {{print $1,$2,$3,$4,$5,$6,$7,$11}}' \\
			| sed 's/Coverage/Coverage\tLoca/' > {outputDir}/amplicon/{sample}.amp.txt
		rm -rf {outputDir}/amplicon/temp
	""".format(sample=sample, outputDir=outputDir, bedtools=bedtools, bamdst=bamdst, bedFile=bedFile)
	os.system(cmd)

# 单肿瘤样本分析
## msisensor
def msisensor(sample, outputDir, panel):
	msisensor = "msisensor"
	msiList = "hg19.msi.list"
	bedFile = "hg19.msisensor.bed"
	if panel == "2":
		bedFile = "hg19.msisensor.test.bed"
	cmd = """
		mkdir {outputDir}/msi
		mkdir {outputDir}/msi/temp
		{msisensor} msi -d {msiList} -t {outputDir}/bam/{sample}.bam \\
			-e {bedFile} -o {outputDir}/msi/temp/{sample}
		cp {outputDir}/msi/temp/{sample} {outputDir}/msi/{sample}.msisensor.txt
		rm -rf {outputDir}/msi/temp
	""".format(sample=sample, outputDir=outputDir, msisensor=msisensor, msiList=msiList, bedFile=bedFile)
	os.system(cmd)

## msisensor2
def msisensor2(sample, outputDir):
	msisensor2 = "msisensor2"
	msiModel = "models_hg19_GRCh37"
	cmd = """
		mkdir {outputDir}/msi/temp
		{msisensor2} msi -M {msiModel} -t {outputDir}/bam/{sample}.bam \\
			-o {outputDir}/msi/temp/{sample}
		cp {outputDir}/msi/temp/{sample} {outputDir}/msi/{sample}.msisensor2.txt
		rm -rf {outputDir}/msi/temp
	""".format(sample=sample, outputDir=outputDir, msisensor2=msisensor2, msiModel=msiModel)
	os.system(cmd)

# visualMSI
def visualMSI(sample, outputDir, panel):
	visualMSI = "visualmsi"
	reference = "hg19.fa"
	bedFile = "hg19.visualmsi.bed"
	if panel == "2":
		bedFile = "hg19.visualmsi.test.bed"
	cmd = """
		mkdir {outputDir}/msi/temp
		{visualMSI} -i {outputDir}/bam/{sample}.bam \\
			-r {reference} -t {bedFile} -j {outputDir}/msi/temp/{sample}.json \\
			-h {outputDir}/msi/temp/{sample}.html > {outputDir}/msi/{sample}.visualmsi.txt
		rm -rf {outputDir}/msi/temp
	""".format(visualMSI=visualMSI, sample=sample, outputDir=outputDir, reference=reference, bedFile=bedFile)
	os.system(cmd)

def tumorOnly(sample, outputDir, panel):
	msisensor(sample, outputDir, panel)
	msisensor2(sample, outputDir)
	visualMSI(sample, outputDir, panel)

# 配对样本分析
## msisensor
def msisensor_pair(sample1, sample2, outputDir, panel):
	msisensor = "msisensor"
	msiList = "hg19.msi.list"
	bedFile = "hg19.msisensor.bed"
	if panel == "2":
		bedFile = "hg19.msisensor.test.bed"
	cmd = """
		if [ ! -d "{outputDir}/msi" ]; then
			mkdir {outputDir}/msi
		fi
		mkdir {outputDir}/msi/temp
		{msisensor} msi -d {msiList} -t {outputDir}/bam/{sample1}.bam \\
			-n {outputDir}/bam/{sample2}.bam \\
			-e {bedFile} -o {outputDir}/msi/temp/{sample1}_{sample2}
		cp {outputDir}/msi/temp/{sample1}_{sample2} {outputDir}/msi/{sample1}_{sample2}.msisensor.txt
		rm -rf {outputDir}/msi/temp
	""".format(sample1=sample1, sample2=sample2, outputDir=outputDir, msisensor=msisensor, msiList=msiList, bedFile=bedFile)
	os.system(cmd)

## visualMSI
def visualMSI_pair(sample1, sample2, outputDir, panel):
	visualMSI = "visualmsi"
	reference = "hg19.fa"
	bedFile = "hg19.visualmsi.bed"
	if panel == "2":
		bedFile = "hg19.visualmsi.test.bed"
	cmd = """
		mkdir {outputDir}/msi/temp
		{visualMSI} -i {outputDir}/bam/{sample1}.bam -n {outputDir}/bam/{sample2}.bam \\
			-r {reference} -t {bedFile} -j {outputDir}/msi/temp/{sample1}_{sample2}.json \\
			-h {outputDir}/msi/temp/{sample1}_{sample2}.html > {outputDir}/msi/{sample1}_{sample2}.visualmsi.txt
		rm -rf {outputDir}/msi/temp
	""".format(visualMSI=visualMSI, sample1=sample1, sample2=sample2, outputDir=outputDir, reference=reference, bedFile=bedFile)
	os.system(cmd)

## mantis
def mantis(sample1, sample2, outputDir, panel):
	mantis = "mantis.py"
	bedFile = "hg19.mantis.bed"
	if panel == "2":
		bedFile = "hg19.mantis.test.bed"
	reference = "hg19.fa"
	cmd = """
		mkdir {outputDir}/msi/temp
		python {mantis} --bedfile {bedFile} --genome {reference} \\
			-n {outputDir}/bam/{sample2}.bam \\
			-t {outputDir}/bam/{sample1}.bam \\
			-o {outputDir}/msi/temp/{sample1}_{sample2}.txt --threads 16
		cp {outputDir}/msi/temp/{sample1}_{sample2}.txt {outputDir}/msi/{sample1}_{sample2}.mantis.txt
		rm -rf {outputDir}/msi/temp
	""".format(sample1=sample1, sample2=sample2, outputDir=outputDir, mantis=mantis, bedFile=bedFile, reference=reference)
	os.system(cmd)

def pairSample(tumorSample, normalSample, outputDir, panel):
	msisensor_pair(tumorSample, normalSample, outputDir, panel)
	visualMSI_pair(tumorSample, normalSample, outputDir, panel)
	mantis(tumorSample, normalSample, outputDir, panel)


# 主流程
def main(sample1, sample2, rawdataDir, outputDir, panel, QC):
	cmd = """
		if [ ! -d "{outputDir}/results" ]; then
			mkdir -p {outputDir}/results
		fi
	""".format(outputDir=outputDir)
	os.system(cmd)

	# 单样本模式
	if sample2 == "none":
		if QC == "yes":
			print "[filter] begin"
			filter(sample1, rawdataDir, outputDir)
		print "[mapping] begin"
		mapping(sample1, rawdataDir, outputDir, QC)
		print "[amplicon] begin"
		amplicon(sample1, outputDir, panel)
		print "[msi analysis] begin"
		tumorOnly(sample1, outputDir, panel)

		results = open(outputDir + "/results/" + sample1 + ".txt", "w")
		msisensorResults = open(outputDir + "/msi/" + sample1 + ".msisensor.txt", "r")
		msisensor2Results = open(outputDir + "/msi/" + sample1 + ".msisensor2.txt", "r")
		visualMSIResults = open(outputDir + "/msi/" + sample1 + ".visualmsi.txt", "r")

		results.write(sample1 + "\t及格位点总数\t不稳定位点数\t百分比\t结论\n")
		for line in msisensorResults:
			if line.startswith("Total"):
				continue
			else:
				lines = line.split("\t")
				totalNumber = lines[0]
				unstable = lines[1]

				if totalNumber == "0":
					results.write("msisensor\t0\t0\t0\t数据质量不足判断\n")
				else:
					percen = str(float(unstable) / float(totalNumber))
					results_o = "\t".join([str(totalNumber), str(unstable), percen])
					if unstable == "0":
						results.write("msisensor\t" + results_o + "\tMSS\n")
					elif unstable == "1":
						results.write("msisensor\t" + results_o + "\tMSI-L\n")
					else:
						results.write("msisensor\t" + results_o + "\tMSI-H\n")
		msisensorResults.close()

		for line in msisensor2Results:
			if line.startswith("Total"):
				continue
			else:
				lines = line.split("\t")
				totalNumber = lines[0]
				unstable = lines[1]

				if totalNumber == "0":
					results.write("msisensor2\t0\t0\t0\t数据质量不足判断\n")
				else:
					percen = str(float(unstable) / float(totalNumber))
					results_o = "\t".join([str(totalNumber), str(unstable), percen])
					if unstable == "0":
						results.write("msisensor2\t" + results_o + "\tMSS\n")
					elif unstable == "1":
						results.write("msisensor2\t" + results_o + "\tMSI-L\n")
					else:
						results.write("msisensor2\t" + results_o + "\tMSI-H\n")
		msisensor2Results.close()
		
		### ！！！！！！设定阈值！！！！！！
		entropyCutOff = 4.2
		
		unstableVisual = 0
		locaVisual = 0
		for line in visualMSIResults:
			if line.startswith("entropy"):
				entropy = line.split("entropy of tumor data: ")[1].split("\n")[0]
				if float(entropy) >= entropyCutOff:
					unstableVisual += 1
			if line.startswith("quality control"):
				if "passed" in line:
					locaVisual += 1
		if locaVisual == 0:
			results.write("visualMSI\t0\t0\t0\t数据质量不足判断\n")
		else:
			percenVisual = str(float(unstableVisual) / float(locaVisual))
			results_o = "\t".join([str(locaVisual), str(unstableVisual), percenVisual])
			if unstableVisual == 0:
				results.write("visualMSI\t" + results_o + "\tMSS\n")
			elif unstableVisual == 1:
				results.write("visualMSI\t" + results_o + "\tMSI-L\n")
			else:
				results.write("visualMSI\t" + results_o + "\tMSI-H\n")
		visualMSIResults.close()
		results.close()

	# 配对样本模式
	else:
		if QC == "yes":
			print "[filter] begin"
			filter(sample1, rawdataDir, outputDir)
			filter(sample2, rawdataDir, outputDir)
		print "[mapping] begin"
		mapping(sample1, rawdataDir, outputDir, QC)
		mapping(sample2, rawdataDir, outputDir, QC)
		print "[amplicon] begin"
		amplicon(sample1, outputDir, panel)
		amplicon(sample2, outputDir, panel)
		print "[msi analysis] begin"
		pairSample(sample1, sample2, outputDir, panel)

		results = open(outputDir + "/results/" + sample1 + "_" + sample2 + ".txt", "w")
		msisensorResults = open(outputDir + "/msi/" + sample1 + "_" + sample2 + ".msisensor.txt", "r")
		visualMSIResults = open(outputDir + "/msi/" + sample1 + "_" + sample2 + ".visualmsi.txt", "r")
		mantisResults = open(outputDir + "/msi/" + sample1 + "_" + sample2 + ".mantis.txt", "r")

		results.write(sample1 + "_" + sample2 + "\t及格位点总数\t不稳定位点数\t百分比\t结论\n")
		for line in msisensorResults:
			if line.startswith("Total"):
				continue
			else:
				lines = line.split("\t")
				totalNumber = lines[0]
				unstable = lines[1]

				if totalNumber == "0":
					results.write("msisensor\t0\t0\t0\t数据质量不足判断\n")
				else:
					percen = str(float(unstable) / float(totalNumber))
					results_o = "\t".join([str(totalNumber), str(unstable), percen])
					if unstable == "0":
						results.write("msisensor\t" + results_o + "\tMSS\n")
					elif unstable == "1":
						results.write("msisensor\t" + results_o + "\tMSI-L\n")
					else:
						results.write("msisensor\t" + results_o + "\tMSI-H\n")
		msisensorResults.close()

		### ！！！！！！设定阈值！！！！！！
		EmdCutOff = 1.0
		
		unstableVisual = 0
		locaVisual = 0
		for line in visualMSIResults:
			if line.startswith("the earth mover"):
				emd = line.split("mover's distance (EMD): ")[1].split("\n")[0]
				if float(emd) >= EmdCutOff:
					unstableVisual += 1
			if line.startswith("quality control"):
				if "passed" in line:
					locaVisual += 1
		if locaVisual == 0:
			results.write("visualMSI\t0\t0\t0\t数据质量不足判断\n")
		else:
			percenVisual = str(float(unstableVisual) / float(locaVisual))
			results_o = "\t".join([str(locaVisual), str(unstableVisual), percenVisual])
			if unstableVisual == 0:
				results.write("visualMSI\t" + results_o + "\tMSS\n")
			elif unstableVisual == 1:
				results.write("visualMSI\t" + results_o + "\tMSI-L\n")
			else:
				results.write("visualMSI\t" + results_o + "\tMSI-H\n")
		visualMSIResults.close()


		# !!!!!!! mantis 阈值 !!!!!!!!
		differenceMantisCutOff = 0.40

		unstableMantis = 0
		locaMantis = 0
		for line in mantisResults:
			if line.startswith("Locus"):
				continue
			elif line.startswith("Average"):
				continue
			else:
				lines = line.split("\t")
				normalReads = int(lines[1])
				tumorReads = int(lines[2])
				diff = float(lines[3])

				if normalReads >= 20 and tumorReads >= 20:
					locaMantis += 1

					if diff >= differenceMantisCutOff:
						unstableMantis += 1
		if locaMantis == 0:
			results.write("mantis\t0\t0\t0\t数据质量不足判断\n")
		else:
			percenMantis = str(float(unstableMantis) / float(locaMantis))
			results_o = "\t".join([str(locaMantis), str(unstableMantis), percenMantis])
			if unstableMantis == 0:
				results.write("mantis\t" + results_o + "\tMSS\n")
			elif unstableMantis == 1:
				results.write("mantis\t" + results_o + "\tMSI-L\n")
			else:
				results.write("mantis\t" + results_o + "\tMSI-H\n")
		mantisResults.close()
		results.close()


if __name__ == "__main__":
	parser = argparse.ArgumentParser(
		description="MSI multiAnalysor PZW",
		prog="MSIpro.py",
		usage="python MSIpro.py [-h] -t <tumor sample id> [-n <normal sample id>] -r <rawdata directory> -o <output directory>",
		formatter_class=argparse.RawTextHelpFormatter
	)
	parser.add_argument("-v", "--version", action="version",
		version="Version 0.2 20200317")
	parser.add_argument("-t", "--tumor", type=str,
		help="Input the tumor sample id")
	parser.add_argument("-n", "--normal", type=str,
		help="Input the normal sample id, option", default="none")
	parser.add_argument("-r", "--rawdata", type=str,
		help="rawdata directory path")
	parser.add_argument("-o", "--output", type=str,
		help="output directory path")
	parser.add_argument("-p", "--panel", type=str,
		help="specify a panel, default = '1' refer to '8 loca panel' and select '2' for '9 loca panel'", default="1")
	parser.add_argument("-q", "--QC", type=str,
		help="do QC, default is 'yes', use '-q no' to skip QC", default="yes")
	if len(sys.argv[1:]) == 0:
		parser.print_help()
		parser.exit()
	args = parser.parse_args()
	main(sample1=args.tumor, sample2=args.normal, rawdataDir=args.rawdata, outputDir=args.output, panel=args.panel, QC=args.QC)