# coding=utf-8
# pzw
# 20200103

import os
import sys
import time
import function.WordWriter3 as ww
import argparse
import json

# 使用软件列表
def softwareDict(software):
	softwareDictionary = {}
	softwareDictionary["bamdst"] = "bamdst"
	softwareDictionary["bedtools"] = "bedtools"
	softwareDictionary["fastp"] = "fastp"
	softwareDictionary["bwa"] = "bwa"
	softwareDictionary["samtools"] = "samtools"
	softwareDictionary["gatk"] = "gatk"

	return softwareDictionary[software]

# 使用数据库列表
def databaseDict(database):
	databaseDictionary = {}
	databaseDictionary["hg19"] = "hg19.fa"
	databaseDictionary["1000G_indels"] = "1000G_phase1.indels.hg19.sites.vcf"
	databaseDictionary["Mills"] = "Mills_and_1000G_gold_standard.indels.hg19.sites.vcf"

	return databaseDictionary[database]

# 获得当前脚本路径
def getAbsPath():
	now = os.path.abspath(os.path.dirname(sys.argv[0]))
	return now

# 获取当前时间
def getTime():
	localtime = time.asctime(time.localtime(time.time()))
	return localtime


# 数据质控
def fastp(rawdataDir, outputDir, sample, thread):
	fastp = softwareDict("fastp")
	cmd = """
		mkdir -p {outputDir}/cleandata
		mkdir {outputDir}/QC
		{fastp} -i {rawdataDir}/{sample}_R1_001.fastq.gz -I {rawdataDir}/{sample}_R2_001.fastq.gz \\
			-o {outputDir}/cleandata/{sample}_1.fq.gz -O {outputDir}/cleandata/{sample}_2.fq.gz \\
			-w {thread} -h {outputDir}/QC/{sample}.html -j {outputDir}/QC/{sample}.json
	""".format(outputDir=outputDir, fastp=fastp, rawdataDir=rawdataDir, sample=sample, thread=thread)

	os.system(cmd)

# 比对
def bwa(outputDir, sample, thread):
	bwa = softwareDict("bwa")
	hg19 = databaseDict("hg19")
	samtools = softwareDict("samtools")
	cmd = """
		mkdir {outputDir}/bam
		{bwa} mem -R "@RG\\tLB:CVD_GOUT\\tPL:illumina\\tPU:CVD_GOUT\\tID:{sample}\\tSM:{sample}" \\
			-t {thread} {hg19} {outputDir}/cleandata/{sample}_1.fq.gz {outputDir}/cleandata/{sample}_2.fq.gz \\
			| {samtools} view -bSh - | {samtools} sort - -@ {thread} -o {outputDir}/bam/{sample}.bam
	""".format(outputDir=outputDir, bwa=bwa, thread=thread, hg19=hg19, sample=sample, samtools=samtools)

	os.system(cmd)

# bam文件校对
def recalibrator(outputDir, sample):
	gatk = softwareDict("gatk")
	samtools = softwareDict("samtools")
	hg19 = databaseDict("hg19")
	indel1 = databaseDict("1000G_indels")
	indel2 = databaseDict("Mills")

	cmd = """
		mkdir {outputDir}/temp
		{gatk} MarkDuplicates -I {outputDir}/bam/{sample}.bam -O {outputDir}/temp/{sample}.md.bam -M {outputDir}/temp/{sample}.md.txt
		{gatk} BuildBamIndex -I {outputDir}/temp/{sample}.md.bam
		{gatk} BaseRecalibrator -I {outputDir}/temp/{sample}.md.bam \\
			-R {hg19} \\
			--known-sites {indel1} \\
			--known-sites {indel2} \\
			-O {outputDir}/temp/{sample}.recal.table
		{gatk} ApplyBQSR -R {hg19} \\
			-I {outputDir}/temp/{sample}.md.bam \\
			--bqsr {outputDir}/temp/{sample}.recal.table \\
			-O {outputDir}/bam/{sample}.final.bam
		{samtools} index {outputDir}/bam/{sample}.final.bam
		rm -rf {outputDir}/temp
	""".format(outputDir=outputDir, gatk=gatk, sample=sample, hg19=hg19, indel1=indel1, indel2=indel2, samtools=samtools)

	os.system(cmd)

# 变异检测
def HaplotypeCaller(outputDir, sample, thread):
	gatk = softwareDict("gatk")
	hg19 = databaseDict("hg19")

	cmd = """
		mkdir {outputDir}/vcf
		{gatk} HaplotypeCaller \\
			-R {hg19} \\
			-I {outputDir}/bam/{sample}.final.bam \\
			--native-pair-hmm-threads {thread} \\
			-O {outputDir}/vcf/{sample}.vcf
	""".format(outputDir=outputDir, gatk=gatk, hg19=hg19, sample=sample, thread=thread)

	os.system(cmd)


# vcf信息获取 ## database = "CVDRiskDB" or "GOUTRiskDB"
def vcfStat(outputDir, sample, database):
	folder = os.path.exists(outputDir + "/insert")
	if not folder:
		os.makedirs(outputDir + "/insert")
	hotspot = open(getAbsPath() + "/data/" + database + ".txt", "r")
	outputFile = open(outputDir + "/insert/" + sample + "." + database +".txt", "w")
	resultsDictFile = open(outputDir + "/insert/" + sample + "." + database + ".results.txt", "w")
	for spot in hotspot:
		if spot.startswith("#"):
			continue
		else:
			spots = spot.split("\t")
			chrom_qurey = spots[2]
			pos_query = spots[3]
			ref_query = spots[4]
			alt_query = spots[5]
			if len(alt_query) >= 3:
				alt_query = "L"
			if len(ref_query) >= 3:
				alt_query = "D"

			type_query = "野生型"
			baseType = ref_query + "/" + ref_query

			vcfFile = open(outputDir + "/vcf/" + sample + ".vcf", "r")
			for line in vcfFile:
				if line.startswith("#"):
					continue
				else:
					lines = line.split("\t")
					chromosome = lines[0]
					pos = lines[1]
					infos = lines[9].split(":")
					genotype = infos[0]
					ADs = infos[1].split(",")
					if len(ADs) != 2:
						continue
					ref_count = ADs[0]
					alt_count = ADs[1]

					if chromosome == chrom_qurey:
						if pos == pos_query:
							# 设定阳性变异阈值
							if (int(ref_count) + int(alt_count)) >= 20:
								if genotype == "1/1":
									type_query = "纯合突变"
									baseType = alt_query + "/" + alt_query
								elif genotype == "0/1":
									type_query = "杂合突变"
									baseType = ref_query + "/" + alt_query
								else:
									continue
			vcfFile.close()
			outputList = [baseType, type_query]
			outputFile.write("\t".join(outputList) + "\n")
			resultList = [chrom_qurey + ":" + pos_query, "".join(baseType.split("/"))]
			resultsDictFile.write("\t".join(resultList) + "\n")
	outputFile.close()
	resultsDictFile.close()
	hotspot.close()

# 风险计算模型
def OrBaseRiskCalculation(resultsDict, database):
	db = open(database, "r")
	dbDict = {}
	for line in db:
		if line.startswith("#"):
			continue
		else:
			lines = line.split("\n")[0].split("\t")
			chrom = lines[2]
			pos = lines[3]
			alt = lines[5]
			OR = float(lines[10])
			frequence = float(lines[9])
			poskey = chrom + ":" + pos
			w = (frequence ** 2) * (OR ** 2) + 2 * frequence * (1 - frequence) * OR + ((1 - frequence) ** 2)
			dbDict[poskey] = [alt, OR, w]
	db.close()
	score_total = 1
	score_normal = 1
	for i in resultsDict.keys():
		n = resultsDict[i].count(dbDict[i][0])
		score = (float(dbDict[i][1]) ** n) / dbDict[i][2]
		score_normal = score_normal * (1 / dbDict[i][2])
		score_total = score_total * score

	risk_power = score_total / score_normal

	return score_total

def main(samplefile, rawdataDir, outputDir, sample, thread):
	sampleDicts = sampleinfo(samplefile)
	CVD_Low_Risk_CutOff = 0.10
	CVD_High_Risk_CutOff = 0.50
	GOUT_Low_Risk_CutOff = 0.10
	GOUT_High_Risk_CutOff = 0.50
	
	print("[样本 " + sample + " 开始]" + getTime())
	print("[过滤开始] " + getTime())
	fastp(rawdataDir, outputDir, sample, thread)
	print("[过滤结束] " + getTime())
	print("[比对开始] " + getTime())
	bwa(outputDir, sample, thread)
	print("[比对结束] " + getTime())
	print("[校对开始] " + getTime())
	recalibrator(outputDir, sample)
	print("[校对结束] " + getTime())
	print("[变异检测开始] " + getTime())
	HaplotypeCaller(outputDir, sample, thread)
	print("[变异检测结束] " + getTime())

	## 报告所需内容提取
	print("[提取内容] " + getTime())
	vcfStat(outputDir, sample, "CVDRiskDB")
	vcfStat(outputDir, sample, "GOUTRiskDB")
	print("[提取内容完成] " + getTime())

	## 风险值计算
	print("[风险值计算] " + getTime())
	CVD_dict = {}
	CVD_dict_File = open(outputDir + "/insert/" + sample + ".CVDRiskDB.results.txt", "r")
	for cvd in CVD_dict_File:
		cvd_k = cvd.split("\n")[0].split("\t")[0]
		cvd_v = cvd.split("\n")[0].split("\t")[1]
		CVD_dict[cvd_k] = cvd_v
	CVD_dict_File.close()

	GOUT_dict = {}
	GOUT_dict_File = open(outputDir + "/insert/" + sample + ".GOUTRiskDB.results.txt", "r")
	for gout in GOUT_dict_File:
		gout_k = gout.split("\n")[0].split("\t")[0]
		gout_v = gout.split("\n")[0].split("\t")[1]
		GOUT_dict[gout_k] = gout_v
	GOUT_dict_File.close()

	CVD_risk = OrBaseRiskCalculation(CVD_dict, getAbsPath() + "/data/CVDRiskDB.txt")
	GOUT_risk = OrBaseRiskCalculation(GOUT_dict, getAbsPath() + "/data/GOUTRiskDB.txt")

	print("CVD风险值为： " + str(CVD_risk))
	print("GOUT风险值为： " + str(GOUT_risk))
	print("[风险值计算完成] " + getTime())

	# 风险阈值判断
	print("[报告填充] " + getTime())
	advice = open(getAbsPath() + "/data/advice.txt", "r")
	adviceDict = {}
	for ad in advice:
		ads = ad.split("\t")
		adviceDict[ads[0]] = ads[1]
	advice.close()

	if CVD_risk <= CVD_Low_Risk_CutOff:
		print("CVD低于平均风险")
		CVD_level = "低于平均水平"
		CVD_image_input = getAbsPath() + "/template/below.jpg"
		CVD_advice = adviceDict["CVD_Low_Risk"]
	elif CVD_risk <= CVD_High_Risk_CutOff:
		print("CVD平均风险")
		CVD_level = "平均水平"
		CVD_image_input = getAbsPath() + "/template/average.jpg"
		CVD_advice = adviceDict["CVD_Average_Risk"]
	else:
		print("CVD高于平均风险")
		CVD_level = "高于平均水平"
		CVD_image_input = getAbsPath() + "/template/above.jpg"
		CVD_advice = adviceDict["CVD_High_Risk"]


	if GOUT_risk <= GOUT_Low_Risk_CutOff:
		print("GOUT低于平均风险")
		GOUT_level = "低于平均水平"
		GOUT_image_input = getAbsPath() + "/template/below.jpg"
		GOUT_advice = adviceDict["GOUT_Low_Risk"]
	elif GOUT_risk <= GOUT_High_Risk_CutOff:
		print("GOUT平均风险")
		GOUT_level = "平均水平"
		GOUT_image_input = getAbsPath() + "/template/average.jpg"
		GOUT_advice = adviceDict["GOUT_Average_Risk"]
	else:
		print("GOUT高于平均风险")
		GOUT_level = "高于平均水平"
		GOUT_image_input = getAbsPath() + "/template/above.jpg"
		GOUT_advice = adviceDict["GOUT_High_Risk"]


	# 报告填写
	fillDict = {}
	fillDict["#[TBS-NAME]#"] = sampleDicts[sample]["name"]
	fillDict["#[TBS-AGE]#"] = sampleDicts[sample]["age"]
	fillDict["#[TBS-SAMPLEID]#"] = sample
	fillDict["#[TBS-FROM]#"] = sampleDicts[sample]["company"]
	fillDict["#[TBS-DOCTOR]#"] = sampleDicts[sample]["doctor"]
	fillDict["#[TBS-DIAGNOISE]#"] = sampleDicts[sample]["diagnosis"]
	fillDict["#[TBS-ILLHISTORY]#"] = sampleDicts[sample]["history"]
	fillDict["#[TBS-DRUGHISTORY]#"] = sampleDicts[sample]["drug"]
	fillDict["#[TBS-GENDER]#"] = sampleDicts[sample]["gender"]
	fillDict["#[TBS-HOSPITALID]#"] = sampleDicts[sample]["hospital_id"]
	fillDict["#[TBS-COLLECTDATE]#"] = sampleDicts[sample]["sampleTime"]
	fillDict["#[TBS-HOSPITALPART]#"] = sampleDicts[sample]["hospital_part"]
	fillDict["#[TBS-PHONE]#"] = sampleDicts[sample]["phone"]
	fillDict["#[HEADER-date]#"] = time.strftime("%Y年%m月%d日", time.localtime())
	fillDict["#[HEADER-name]#"] = sampleDicts[sample]["name"]

	fillDict["#[CVD-risk]#"] = str("%.2f" % CVD_risk)
	fillDict["#[CVD-risklevel]#"] = CVD_level
	fillDict["#[IMAGE-CVD-(40,10)]#"] = CVD_image_input
	fillDict["#[CVD-advice]#"] = CVD_advice
	fillDict["#[TABLE-CVD]#"] = outputDir + "/insert/" + sample + ".CVDRiskDB.txt"

	fillDict["#[GOUT-risk]#"] = str("%.2f" % GOUT_risk)
	fillDict["#[GOUT-risklevel]#"] = GOUT_level
	fillDict["#[IMAGE-GOUT-(40,10)]#"] = GOUT_image_input
	fillDict["#[GOUT-advice]#"] = GOUT_advice
	fillDict["#[TABLE-GOUT]#"] = outputDir + "/insert/" + sample + ".GOUTRiskDB.txt"

	folder = os.path.exists(outputDir + "/report")
	if not folder:
		os.makedirs(outputDir + "/report")

	outputdocx = outputDir + "/report/" + sample + "-CVD_GOUT-" + fillDict["#[TBS-NAME]#"] + "-" + fillDict["#[HEADER-date]#"] + ".docx"
	ww.WordWriter(getAbsPath() + "/template/CVD_20191230.docx", outputdocx, fillDict)
	print("[报告填充完成] " + getTime())
	print("[样本" + sample + "完成] " + getTime())


if __name__ == "__main__":
	parser = argparse.ArgumentParser(
		description="CVD & GOUT Risk Calculator and Reporter",
		prog="CVDRiskAnalysis.py",
		usage="python3 CVDRiskAnalysis.py [-h] -s <sampleinfo> -r <rawdataDir> -o <outputDir> -i <sampleID> -t <threads>",
		formatter_class=argparse.RawTextHelpFormatter
	)
	parser.add_argument("-v", "--version", action="version",
		version="Version 0.2 20200103")
	parser.add_argument("-s", "--sampleinfo", type=str,
		help="Input the sampleinfo file", default=getAbsPath() + "/data/CVD.txt")
	parser.add_argument("-r", "--rawdata", type=str,
		help="Input the rawdata directory")
	parser.add_argument("-o", "--output", type=str,
		help="specify a output directory")
	parser.add_argument("-i", "--sample", type=str,
		help="give me a sample ID", default="CVD")
	parser.add_argument("-t", "--threads", type=str,
		help="threads to run, default=4", default="4")

	if len(sys.argv[1:]) == 0:
		parser.print_help()
		parser.exit()
	args = parser.parse_args()
	main(samplefile=args.sampleinfo, rawdataDir=args.rawdata, outputDir=args.output, sample=args.sample, thread=args.threads)