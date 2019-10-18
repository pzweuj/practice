# pzw
# 20190929
# CLOCK ANALYSIS PIPELINE

import os
import json
import pandas as pd
import numpy as np
import time
import argparse
import sys

def QC(sample, rawdataDir, outputDir):
	cmd = """
		if [ ! -d "{outputDir}" ]; then
			mkdir {outputDir}
		fi
		mkdir {outputDir}/cleandata
		mkdir {outputDir}/QC
		fastp -i {rawdataDir}/{sample}_R1_001.fastq.gz \\
			-I {rawdataDir}/{sample}_R2_001.fastq.gz \\
			-o {outputDir}/cleandata/{sample}_1.fq.gz \\
			-O {outputDir}/cleandata/{sample}_2.fq.gz \\
			-w 8 -j {outputDir}/QC/{sample}.json \\
			-h {outputDir}/QC/{sample}.html
	""".format(sample=sample, rawdataDir=rawdataDir, outputDir=outputDir)
	os.system(cmd)

def mapping(sample, outputDir):
	cmd = """
		mkdir {outputDir}/bam
		bwa mem -t 8 \\
			-R "@RG\\tID:{sample}\\tSM:{sample}\\tLB:CLOCKPanel\\tPL:illumina" \\
			/home/zhaowen/workspace/database/human/hg19.fa \\
			{outputDir}/cleandata/{sample}_1.fq.gz \\
			{outputDir}/cleandata/{sample}_2.fq.gz \\
			| samtools view -bSh - \\
			| samtools sort -@ 8 - -o {outputDir}/bam/{sample}.bam

	""".format(sample=sample, outputDir=outputDir)
	os.system(cmd)

def targetBam(sample, outputDir):
	cmd = """
		mkdir {outputDir}/bed
		bedtools genomecov -ibam {outputDir}/bam/{sample}.bam -bga \\
			| bedtools intersect -a - \\
			-b /home/zhaowen/workspace/database/human/clock.bed -wb \\
			> {outputDir}/bed/{sample}.txt
		bedtools intersect -a {outputDir}/bam/{sample}.bam \\
			-b /home/zhaowen/workspace/database/human/clock.bed \\
			> {outputDir}/bam/{sample}.target.bam
		samtools view -F 2052 {outputDir}/bam/{sample}.target.bam \\
			| wc -l > {outputDir}/bed/{sample}.read.txt
	""".format(sample=sample, outputDir=outputDir)
	os.system(cmd)

def VariantCall(sample, outputDir):
	cmd = """
		mkdir {outputDir}/vcf
		samtools index {outputDir}/bam/{sample}.target.bam
		gatk HaplotypeCaller -R /home/zhaowen/workspace/database/human/hg19.fa \\
			-I {outputDir}/bam/{sample}.target.bam \\
			-O {outputDir}/vcf/{sample}.vcf \\
			--native-pair-hmm-threads 8
	""".format(sample=sample, outputDir=outputDir)
	os.system(cmd)

def anno(sample, outputDir):
	cmd = """
		mkdir {outputDir}/annotation
		/software/annovar2018/convert2annovar.pl -format vcf4 {outputDir}/vcf/{sample}.vcf \\
			-includeinfo > {outputDir}/annotation/{sample}.avinput
		/software/annovar2018/table_annovar.pl {outputDir}/annotation/{sample}.avinput \\
			/software/annovar2018/humandb -buildver hg19 \\
			-out {outputDir}/annotation/{sample} -remove \\
			-protocol refGene,cytoBand,1000g2015aug_all,1000g2015aug_eas,avsnp150,clinvar_20190305,cosmic70,dbnsfp35a \\
			-operation g,r,f,f,f,f,f,f \\
			-nastring . -thread 8 -otherinfo
		rm {outputDir}/annotation/{sample}.avinput
	""".format(sample=sample, outputDir=outputDir)
	os.system(cmd)

def sv_detect(sample, outputDir):
	cmd = """
		samtools index {outputDir}/bam/{sample}.bam
		samtools view -bh -F 1294 {outputDir}/bam/{sample}.bam \\
			| samtools sort -@ 8 -o {outputDir}/bam/{sample}.discordants.bam
		samtools view -h {outputDir}/bam/{sample}.bam \\
			| /home/zhaowen/workspace/software/lumpy-sv/scripts/extractSplitReads_BwaMem \\
			-i stdin | samtools view -bSh - \\
			| samtools sort -@ 8 - -o {outputDir}/bam/{sample}.splitters.bam
		samtools index {outputDir}/bam/{sample}.discordants.bam
		samtools index {outputDir}/bam/{sample}.splitters.bam
		lumpyexpress -B {outputDir}/bam/{sample}.bam \\
			-S {outputDir}/bam/{sample}.splitters.bam \\
			-D {outputDir}/bam/{sample}.discordants.bam \\
			-o {outputDir}/vcf/{sample}.sv.vcf
		bedtools intersect -a {outputDir}/vcf/{sample}.sv.vcf \\
			-b /home/zhaowen/workspace/database/human/clock.bed -header \\
			> {outputDir}/vcf/{sample}.sv.target.vcf
		svtyper-sso -i {outputDir}/vcf/{sample}.sv.target.vcf \\
			-B {outputDir}/bam/{sample}.bam \\
			--cores 8 -o {outputDir}/vcf/{sample}.gt.target.vcf
	""".format(sample=sample, outputDir=outputDir)
	os.system(cmd)

def jsonFileEx(sample, outputDir):
	jsonFile = json.load(outputDir + "/QC/" + sample + ".json")
	rawReads = jsonFile["summary"]["before_filtering"]["total_reads"]
	rawBases = jsonFile["summary"]["before_filtering"]["total_bases"]
	gcContent = jsonFile["summary"]["before_filtering"]["gc_content"]
	cleanBasesQ20 = jsonFile["summary"]["after_filtering"]["q20_bases"]
	cleanBasesQ30 = jsonFile["summary"]["after_filtering"]["q30_bases"]
	return [rawReads, rawBases, gcContent, cleanBasesQ20, cleanBasesQ30]

# def ReadsOnTarget(sample, outputDir):
# 	countFile = open(outputDir + "/bed/" + sample + ".read.txt")
# 	for line in countFile:
# 		if line != "":
# 			targetReads = int(line.split("\n")[0])
# 	countFile.close()
# 	return targetReads

def QCreport(sample, outputDir):

	jsonFile = json.load(open(outputDir+"/QC/"+sample+".json", "r"))

	rawReads = jsonFile["summary"]["before_filtering"]["total_reads"]
	rawBases = str(jsonFile["summary"]["before_filtering"]["total_bases"])
	# rawGC = jsonFile["summary"]["before_filtering"]["gc_content"]
	# rawQ20 = jsonFile["summary"]["before_filtering"]["q20_bases"]
	# rawQ30 = jsonFile["summary"]["before_filtering"]["q30_bases"]
	cleanGC = jsonFile["summary"]["after_filtering"]["gc_content"]
	cleanQ20 = str(jsonFile["summary"]["after_filtering"]["q20_bases"])
	cleanQ30 = str(jsonFile["summary"]["after_filtering"]["q30_bases"])
	cleanReads = jsonFile["summary"]["after_filtering"]["total_reads"]
	cleanBases = jsonFile["summary"]["after_filtering"]["total_bases"]

	for line in open(outputDir+"/bed/"+sample+".read.txt", "r"):
		targetReads = int(line)

	OnTargetReadsPer = str("%.2f" % ((float(targetReads) / float(cleanReads)) * 100)) + "%"

	covData = open(outputDir+"/bed/"+sample+".txt", "r")
	OnTargetBases = 0
	bedLength = 0
	covBases = 0
	covBases10 = 0
	covBases30 = 0
	covBases100 = 0
	for i in covData:
		isplit = i.split("\t")
		start = isplit[1]
		end = isplit[2]
		cov = isplit[3]
		bedLength += int(end) - int(start)
		OnTargetBases += (int(cov) * (int(end) - int(start)))

		if int(cov) > 0:
			covBases += int(end) - int(start)

		if int(cov) >= 10:
			covBases10 += int(end) - int(start)

		if int(cov) >= 30:
			covBases30 += int(end) - int(start)

		if int(cov) >= 100:
			covBases100 += int(end) - int(start)


	OnTargetBasesPer = str("%.2f" % ((float(OnTargetBases) / float(cleanBases)) * 100)) + "%"

	averageDepth = float(OnTargetBases) / float(bedLength)

	covPer = float(covBases) / float(bedLength)
	covPer10 = float(covBases10) / float(bedLength)
	covPer30 = float(covBases30) / float(bedLength)
	covPer100 = float(covBases100) / float(bedLength)

	output = [
		sample,
		str(rawReads),
		rawBases,
		str(cleanReads),
		str(cleanBases),
		str(cleanGC),
		cleanQ20,
		cleanQ30,
		str(targetReads),
		OnTargetReadsPer,
		str(OnTargetBases),
		str(OnTargetBasesPer),
		str(averageDepth),
		str(covPer),
		str(covPer10),
		str(covPer30),
		str(covPer100)
	]

	return output

def geneCov(sample, outputDir):
	geneList = ["PER3", "CORT", "CHRNB2", "NPAS2", "NCKAP5", "PER2", "BHLHE40", "PROK2", "CLOCK", "NOCT", "CSNK1A1", "FAM50B", "KCNT1", "PTGDS", "ARNTL", "CRY2", "RBM4B", "BHLHE41", "ARNTL2", "TIMELESS", "PMCH", "CRY1", "CIPC", "PER1", "NR1D1", "RARA", "HCRT", "AANAT", "CRTC1", "CHRNA4 ", "DEPDC5", "PASD1"]
	outputDict = {}
	for i in geneList:
		covData = open(outputDir+"/bed/"+sample+".txt", "r")

		length = 0
		TargetBases = 0
		covLength = 0
		covLength10 = 0
		covLength30 = 0
		covLength100 = 0
		for line in covData:
			if i in line:
				ii = line.split("\n")[0].split("\t")

				chrom = ii[4]
				geneStart = ii[5]
				geneEnd = ii[6]

				cov = ii[3]
				start = ii[1]
				end = ii[2]

				length += int(end) - int(start)

				TargetBases += (int(cov) * (int(end) - int(start)))

				if int(cov) > 0:
					covLength += int(end) - int(start)
				if int(cov) > 10:
					covLength10 += int(end) - int(start)
				if int(cov) > 30:
					covLength30 += int(end) - int(start)
				if int(cov) > 100:
					covLength100 += int(end) - int(start)

		averageDepth = float(TargetBases) / length
		covPer = float(covLength) / length
		covPer10 = float(covLength10) / length	 
		covPer30 = float(covLength30) / length
		covPer100 = float(covLength100) / length

		outputList = [chrom, geneStart, geneEnd, str(length), str(averageDepth), str(covPer), str(covPer10), str(covPer30), str(covPer100)]
		outputDict[i] = outputList
		covData.close()

	return outputDict

def annovarFix(sample, outputDir):
	annoFile = open(outputDir+"/annotation/"+sample+".hg19_multianno.txt", "r")
	annoFileFix = open(outputDir+"/annotation/"+sample+".fix.txt", "w")
	for line in annoFile:
		if line.startswith("Chr"):
			i = 1
			otherInfoTitleFix = []
			while i <= 9:
				otherInfoFix = "otherinfo" + str(i)
				otherInfoTitleFix.append(otherInfoFix)
				i += 1

			otherInfoFixString = "\t".join(otherInfoTitleFix)

			annoFileFix.write(line.split("\n")[0] + "\t" + otherInfoFixString + "\n")
		else:
			annoFileFix.write(line)
	annoFile.close()
	annoFileFix.close()

	pd.set_option("mode.chained_assignment", None)
	df = pd.read_csv(outputDir+"/annotation/"+sample+".fix.txt", sep="\t", header=0)

	df_exon = df[df["Func.refGene"] == "exonic"]

	df_exon.insert(7, "BaseChange", ".")
	df_exon.insert(8, "AAChange", ".")
	df_exon.insert(9, "exon", ".")
	df_exon.insert(10, "AO", ".")
	df_exon.insert(11, "DP", ".")
	df_exon.insert(12, "AF", ".")

	for i in df_exon.index:
		df_exon.loc[i, "AO"] = df_exon.loc[i, "otherinfo9"].split(":")[1].split(",")[1]
		df_exon.loc[i, "DP"] = df_exon.loc[i, "otherinfo9"].split(":")[2]
		df_exon.loc[i, "AF"] = float(df_exon.loc[i, "AO"]) / float(df_exon.loc[i, "DP"])

		if "," in df_exon.loc[i, "AAChange.refGene"]:
			df_exon.loc[i, "BaseChange"] = df_exon.loc[i, "AAChange.refGene"].split(",")[0].split(":")[3]
			df_exon.loc[i, "AAChange"] = df_exon.loc[i, "AAChange.refGene"].split(",")[0].split(":")[4]
			df_exon.loc[i, "exon"] = df_exon.loc[i, "AAChange.refGene"].split(",")[0].split(":")[2]
		else:
			df_exon.loc[i, "BaseChange"] = df_exon.loc[i, "AAChange.refGene"].split(":")[3]
			df_exon.loc[i, "AAChange"] = df_exon.loc[i, "AAChange.refGene"].split(":")[4]
			df_exon.loc[i, "exon"] = df_exon.loc[i, "AAChange.refGene"].split(":")[2]


	df_exon.drop(["GeneDetail.refGene", "CLNALLELEID", "CLNDISDB", "CLNREVSTAT", "Otherinfo", "otherinfo1", "otherinfo2", "otherinfo3", "otherinfo4", "otherinfo5", "otherinfo6", "otherinfo7", "otherinfo8", "otherinfo9"], axis=1, inplace=True)
	df_exon.to_csv(outputDir+"/annotation/"+sample+".anno.txt", sep="\t", header=True, index=False)
	print "Fix " + sample + " done!"

def svInfo(sample, outputDir):
	svFile = open(outputDir+"/vcf/"+sample+".gt.target.vcf", "r")
	outputDict = {}
	i = 0
	for line in svFile:
		if line.startswith("#"):
			continue
		elif "IMPRECISE" in line:
			continue
		else:
			lines = line.split("\t")
			infos = lines[7]

			infosDetails = infos.split(";")
			GQ = lines[9].split(":")[4]
			SU = lines[9].split(":")[1]
			combineInfo = lines[4]
			if "BND" in infosDetails[0]:
				CIPOS95 = "[" + infosDetails[4].split("=")[1] + "]"
				CIEND95 = "[" + infosDetails[5].split("=")[1] + "]"

				if combineInfo[0] == "N":
					if "[" in combineInfo:
						i += 1
						chrom1 = lines[0]
						breakpoint1 = lines[1]
						strand1 = infosDetails[1].split("=")[1][0]
						svType = combineInfo.replace("N", chrom1+":"+breakpoint1)
						chrom2 = combineInfo.split("[")[1].split(":")[0]
						breakpoint2 = combineInfo.split("[")[1].split(":")[1]
						strand2 = infosDetails[1].split("=")[1][1]
						outputDict[str(i)] = [chrom1, breakpoint1, strand1, chrom2, breakpoint2, strand2, svType, GQ, SU]

					elif "]" in combineInfo:
						i += 1
						chrom1 = lines[0]
						breakpoint1 = lines[1]
						strand1 = infosDetails[1].split("=")[1][0]
						svType = combineInfo.replace("N", chrom1+":"+breakpoint1)
						chrom2 = combineInfo.split("]")[1].split(":")[0]
						breakpoint2 = combineInfo.split("]")[1].split(":")[1]
						strand2 = infosDetails[1].split("=")[1][1]
						outputDict[str(i)] = [chrom1, breakpoint1, strand1, chrom2, breakpoint2, strand2, svType, GQ, SU]

					else:
						continue
				elif combineInfo[-1] == "N":
					if "[" in combineInfo:
						i += 1
						chrom2 = lines[0]
						breakpoint2 = lines[1]
						strand2 = infosDetails[1].split("=")[1][1]
						svType = combineInfo.replace("N", chrom2+":"+breakpoint2)
						chrom1 = combineInfo.split("[")[1].split(":")[0]
						breakpoint1 = combineInfo.split("[")[1].split(":")[1]
						strand1 = infosDetails[1].split("=")[1][0]
						outputDict[str(i)] = [chrom1, breakpoint1, strand1, chrom2, breakpoint2, strand2, svType, GQ, SU]

					elif "]" in combineInfo:
						i += 1
						chrom2 = lines[0]
						breakpoint2 = lines[1]
						strand2 = infosDetails[1].split("=")[1][1]
						svType = combineInfo.replace("N", chrom2+":"+breakpoint2)
						chrom1 = combineInfo.split("]")[1].split(":")[0]
						breakpoint1 = combineInfo.split("]")[1].split(":")[1]
						strand1 = infosDetails[1].split("=")[1][0]
						outputDict[str(i)] = [chrom1, breakpoint1, strand1, chrom2, breakpoint2, strand2, svType, GQ, SU]

					else:
						continue
				else:
					continue
	svFile.close()
	return outputDict

def PrintOutResults(sample, outputDir):
	folder = os.path.exists(outputDir+"/results")
	if not folder:
		os.makedirs(outputDir+"/results")

	# QC
	QCFile = open(outputDir+"/results/"+sample+".QC.txt", "w")
	QCStat = QCreport(sample, outputDir)
	qcOutput = [
		"SampleID",
		"rawReads",
		"rawBases",
		"cleanReads",
		"cleanBases",
		"cleanDataGC",
		"cleanData_Bases_Q20",
		"cleanData_Bases_Q30",
		"Target_Reads",
		"Target_Reads_Percent",
		"Target_Bases",
		"Target_Bases_Percent",
		"average_Depth",
		"coverage",
		"coverage_above_10x",
		"coverage_above_30x",
		"coverage_above_100x"
	]
	QCFile.write("\t".join(qcOutput) + "\n")
	QCFile.write("\t".join(QCStat) + "\n")

	QCFile.close()

	# coverage
	coverFile = open(outputDir+"/results/"+sample+".coverage.txt", "w")
	coverFile.write("\t".join(["Gene", "Chr", "Start", "End", "Length", "Average_Depth", "Coverage", "Coverage_10X", "Coverage_30X", "Coverage_100X"]) + "\n")
	covStat = geneCov(sample, outputDir)
	for i in covStat.keys():
		coverFile.write(i + "\t" + "\t".join(covStat[i]) + "\n")
	coverFile.close()


	# annovar
	annovarFix(sample, outputDir)
	cmd = """
		cp {outputDir}/annotation/{sample}.anno.txt {outputDir}/results/{sample}.finalAnno.txt
	""".format(sample=sample, outputDir=outputDir)
	os.system(cmd)

	# SV
	svFile = open(outputDir+"/results/"+sample+".sv.txt", "w")
	svStat = svInfo(sample, outputDir)
	svFile.write("\t".join(["Chromosome1", "Pos1", "Strand1", "Chromosome2", "Pos2", "Strand2", "svType", "GQ", "ReadsCounts"]) + "\n")
	for i in svStat.keys():
		svFile.write("\t".join(svStat[i]) + "\n")
	svFile.close()

def getTime():
	localtime = time.asctime(time.localtime(time.time()))
	return localtime

def main(sample, rawdataDir, outputDir):
	print "***[ "+sample+" start]***", getTime()
	print "[filter start]", getTime()
	QC(sample, rawdataDir, outputDir)
	print "[mapping start]", getTime()
	mapping(sample, outputDir)
	print "[targetBam start]", getTime()
	targetBam(sample, outputDir)
	print "[VariantCall start]", getTime()
	VariantCall(sample, outputDir)
	print "[Annotate start]", getTime()
	anno(sample, outputDir)
	print "[SV Detect start]", getTime()
	sv_detect(sample, outputDir)
	print "[Output Data start]", getTime()
	PrintOutResults(sample, outputDir)
	print "***[ "+sample+" done]***", getTime()

if __name__ == "__main__":
	parser = argparse.ArgumentParser(description="CLOCK Analysis Pipeline by genephar",
		prog="clock_analysis",
		usage="python clock_analysis.py -i <sample> -o <output dir> -r <rawdata dir>")
	group = parser.add_mutually_exclusive_group()
	parser.add_argument("-v", "--version", action="version",
		version="Version 0.1 20190930")
	parser.add_argument("-i", "--input", type=str,
		help="Input the sample ID")
	parser.add_argument("-o", "--output", type=str,
		help="the output directory")
	parser.add_argument("-r", "--rawdata", type=str,
		help="the rawdata directory")
	if len(sys.argv[1:]) == 0:
		parser.print_help()
		parser.exit()
	args = parser.parse_args()
	main(sample=args.input, outputDir=args.output, rawdataDir=args.rawdata)
