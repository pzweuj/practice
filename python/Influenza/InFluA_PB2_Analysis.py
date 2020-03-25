# coding=utf-8
# pzw
# 20203024

import os
import json
from openpyxl import load_workbook
import datetime
import sys
import argparse

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

def mapping(sample, outputDir):
	reference = "PB2.fasta"
	bwa = "bwa"
	cmd = """
		mkdir {outputDir}/bam
		{bwa} mem -t 8 {reference} {outputDir}/cleandata/{sample}_1.fq.gz \\
			{outputDir}/cleandata/{sample}_2.fq.gz \\
			| samtools view -bSh - | samtools sort -@ 8 - -o {outputDir}/bam/{sample}.bam
		samtools index {outputDir}/bam/{sample}.bam
	""".format(bwa=bwa, reference=reference, outputDir=outputDir, sample=sample)
	os.system(cmd)

def calling(sample, outputDir):
	reference = "PB2.fasta"
	bcftools = "bcftools"
	cmd = """
		mkdir {outputDir}/vcf
		{bcftools} mpileup -f {reference} \\
			{outputDir}/bam/{sample}.bam -d 1000 \\
			| {bcftools} call -mv -O v \\
			-o {outputDir}/vcf/{sample}.vcf
	""".format(outputDir=outputDir, bcftools=bcftools, reference=reference, sample=sample)
	os.system(cmd)


def BamStat(sample, outputDir):
	bamdst = "bamdst"
	bedfile = "databases/PB2.bed"
	cmd = """
		mkdir {outputDir}/temp
		{bamdst} -p {bedfile} {outputDir}/bam/{sample}.bam -o {outputDir}/temp
		zcat {outputDir}/temp/region.tsv.gz > {outputDir}/temp/{sample}.txt
	""".format(outputDir=outputDir, bamdst=bamdst, bedfile=bedfile, sample=sample)
	os.system(cmd)
	stat = open(outputDir + "/temp/" + sample + ".txt", "r")
	for line in stat:
		if line.startswith("#"):
			continue
		else:
			ls = line.split("\t")
			avg_depth = ls[3]
			cov = ls[5]
	cmd2 = """
		rm -rf {outputDir}/temp
	""".format(outputDir=outputDir)
	os.system(cmd2)
	return [avg_depth, cov]

def GetPB2Seq(sample, outputDir):
	PB2 = open("PB2.fasta", "r")
	PB2_seq = ""
	for line in PB2:
		if line.startswith(">"):
			continue
		else:
			PB2_seq = PB2_seq + line.split("\n")[0]
	PB2.close()

	vcfFile = open(outputDir + "/vcf/" + sample + ".vcf", "r")
	changeDict = {}
	for vcf in vcfFile:
		if vcf.startswith("#"):
			continue
		else:
			vcfs = vcf.split("\t")
			pos = vcfs[1]
			alt = vcfs[4]
			DP = vcfs[7].split(";")[0].split("=")[1]
			if int(DP) >= 100:
				changeDict[pos] = alt
	vcfFile.close()

	for loc in changeDict.keys():
		location = int(loc)
		PB2_seq = PB2_seq[:location-1] + changeDict[loc] + PB2_seq[location:]

	return PB2_seq

def translate(seq):
	table = { 
		"ATA":"I", "ATC":"I", "ATT":"I", "ATG":"M",
		"ACA":"T", "ACC":"T", "ACG":"T", "ACT":"T",
		"AAC":"N", "AAT":"N", "AAA":"K", "AAG":"K",
		"AGC":"S", "AGT":"S", "AGA":"R", "AGG":"R",
		"CTA":"L", "CTC":"L", "CTG":"L", "CTT":"L",
		"CCA":"P", "CCC":"P", "CCG":"P", "CCT":"P",
		"CAC":"H", "CAT":"H", "CAA":"Q", "CAG":"Q",
		"CGA":"R", "CGC":"R", "CGG":"R", "CGT":"R",
		"GTA":"V", "GTC":"V", "GTG":"V", "GTT":"V",
		"GCA":"A", "GCC":"A", "GCG":"A", "GCT":"A",
		"GAC":"D", "GAT":"D", "GAA":"E", "GAG":"E",
		"GGA":"G", "GGC":"G", "GGG":"G", "GGT":"G",
		"TCA":"S", "TCC":"S", "TCG":"S", "TCT":"S",
		"TTC":"F", "TTT":"F", "TTA":"L", "TTG":"L",
		"TAC":"Y", "TAT":"Y", "TAA":"_", "TAG":"_",
		"TGC":"C", "TGT":"C", "TGA":"_", "TGG":"W",
	}
	protein = ""
	if len(seq) % 3 == 0:
		for i in range(0, len(seq), 3):
			codon = seq[i: i + 3]
			protein += table[codon]
	else:
		print("Warning, sequence can not be translated completely")
		for i in range(0, len(seq), 3):
			if i + 3 < len(seq):
				codon = seq[i: i + 3]
				protein += table[codon]
	return protein

def GetPB2Protein(seq):
	PB2_seq = seq[27:2304]
	PB2_protein = translate(PB2_seq)
	return PB2_protein

def ImportantLoca(protein):
	databasesLocation = "databases/databases.txt"
	databases = open(databasesLocation, "r")
	importList = []
	for line in databases:
		if line.startswith("#"):
			continue
		else:
			ll = line.split("\n")[0].split("\t")
			loca = (ll[1] + ll[0] + ll[2])
			location = int(ll[0])
			pmid = ll[3]
			if protein[location - 1] == ll[2]:
				result = "positive"
			else:
				result = "-"
			output = [loca, result, pmid]
			importList.append(output)
	return importList

def main(sample, rawdataDir, outputDir):
	filter(sample, rawdataDir, outputDir)
	mapping(sample, outputDir)
	calling(sample, outputDir)

	PB2Sequence = GetPB2Seq(sample, outputDir)
	PB2ProteinSequence = GetPB2Protein(PB2Sequence)

	importantLocaChangeList = ImportantLoca(PB2ProteinSequence)

	jsonFile = json.load(open(outputDir + "/QC/" + sample + ".json", "r"))
	before = jsonFile["summary"]["before_filtering"]
	after = jsonFile["summary"]["after_filtering"]
	before_reads = before["total_reads"]
	before_bases = before["total_bases"]
	before_q20 = before["q20_rate"]
	before_q30 = before["q30_rate"]
	before_gc = before["gc_content"]
	after_reads = after["total_reads"]
	after_bases = after["total_bases"]
	after_q20 = after["q20_rate"]
	after_q30 = after["q30_rate"]
	after_gc = after["gc_content"]
	stats = BamStat(sample, outputDir)

	wb = load_workbook("databases/PB2_20200324.xlsx")
	summary = wb.get_sheet_by_name("summary")
	sequence = wb.get_sheet_by_name("Sequence")
	vcf = wb.get_sheet_by_name("Vcf")

	## summary
	summary["A8"] = sample
	dateStr = "".join(str(datetime.date.today()).split("-"))
	summary["A11"] = dateStr
	summary["A19"] = before_reads
	summary["B19"] = before_bases
	summary["C19"] = before_q20
	summary["D19"] = before_q30
	summary["E19"] = before_gc
	summary["F19"] = after_reads
	summary["G19"] = after_bases
	summary["H19"] = after_q20
	summary["I19"] = after_q30
	summary["J19"] = after_gc
	summary["K19"] = stats[0]
	summary["L19"] = stats[1]

	rows = ImportantLoca(PB2ProteinSequence)
	nLen = len(rows)
	s = 24
	b = 0
	while b < nLen:
		cA = "A" + str(s)
		cB = "B" + str(s)
		cC = "C" + str(s)

		summary[cA] = rows[b][0]
		summary[cB] = rows[b][1]
		summary[cC] = rows[b][2]

		b += 1
		s += 1

	## Sequence
	sequence["A1"] = ">" + sample + "_PB2_origin"
	n = 2
	while n <= 25:
		c = "A" + str(n)
		start = (n - 2) * 100
		end = (n - 1) * 100
		sequence[c] = PB2Sequence[start:end]
		n += 1

	sequence["A28"] = ">" + sample + "_PB2_translation"
	m = 29
	while m <= 36:
		cc = "A" + str(m)
		startm = (m - 29) * 100
		endm = (m - 28) * 100
		sequence[cc] = PB2ProteinSequence[startm:endm]
		m += 1

	## Vcf
	vcf["J1"] = sample
	vcfFile = open(outputDir + "/vcf/" + sample + ".vcf", "r")
	for line in vcfFile:
		if line.startswith("#"):
			continue
		else:
			row = line.split("\n")[0].split("\t")
			vcf.append(row)
	vcfFile.close()

	cmd = """
		mkdir {outputDir}/results
	""".format(outputDir=outputDir)
	os.system(cmd)
	outputName = "PB2_" + sample + "_" + dateStr + ".xlsx"
	wb.save(outputDir + "/results/" + outputName)
	print("Task done!")

if __name__ == "__main__":
	parser = argparse.ArgumentParser(
		description="influenza A PB2 PZW",
		prog="InFluA_PB2_Analysis.py",
		usage="python3 InFluA_PB2_Analysis.py [-h] -r <rawdataDir> -o <outputDir> -i <sampleID>",
		formatter_class=argparse.RawTextHelpFormatter
	)
	parser.add_argument("-v", "--version", action="version",
		version="Version 0.2 20200324")
	parser.add_argument("-r", "--rawdata", type=str,
		help="rawdata directory path")
	parser.add_argument("-o", "--output", type=str,
		help="output directory path", default="none")
	parser.add_argument("-i", "--sample", type=str,
		help="sample ID")
	if len(sys.argv[1:]) == 0:
		parser.print_help()
		parser.exit()
	args = parser.parse_args()
	main(sample=args.sample, rawdataDir=args.rawdata, outputDir=args.output)