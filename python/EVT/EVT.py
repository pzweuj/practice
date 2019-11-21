# coding=utf-8
# EVT pipeline
# pzw
# 20191120

import os
import sys
import bamnostic
from collections import Counter
import argparse
import json
import commands
from Bio.Blast import NCBIWWW
import xml.etree.cElementTree as ET


# 原始数据质控与过滤
def filter(sample, rawdata, directory):
	cmd = """
		mkdir -p {directory}/QC
		mkdir -p {directory}/cleandata

		fastp -i {rawdata}/{sample}_R1_001.fastq.gz -I {rawdata}/{sample}_R2_001.fastq.gz \\
			-o {directory}/cleandata/{sample}_1.fq.gz -O {directory}/cleandata/{sample}_2.fq.gz \\
			-w 8 -j {directory}/QC/{sample}.json -h {directory}/QC/{sample}.html
	""".format(sample=sample, rawdata=rawdata, directory=directory)
	os.system(cmd)

# 比对到人参考基因组hg19，并过滤人源reads
def mapToHost(sample, directory):
	cmd = """
		mkdir -p {directory}/host
		bwa mem -t 8 /home/zhaowen/workspace/database/human/hg19.fa \\
			{directory}/cleandata/{sample}_1.fq.gz \\
			{directory}/cleandata/{sample}_2.fq.gz \\
			| samtools view -bSh - | samtools sort -@ 8 - -o {directory}/host/{sample}.bam
		samtools view {directory}/host/{sample}.bam -F 2052 | wc -l > {directory}/host/{sample}.hostCount.txt
		samtools view -bh {directory}/host/{sample}.bam -f 12 -F 256 \\
			| samtools sort -n -@ 8 - -o {directory}/host/{sample}.nohost.bam
		bedtools bamtofastq -i {directory}/host/{sample}.nohost.bam \\
			-fq {directory}/host/{sample}_1.fq \\
			-fq2 {directory}/host/{sample}_2.fq
		gzip {directory}/host/{sample}_1.fq
		gzip {directory}/host/{sample}_2.fq
		rm {directory}/host/{sample}.bam {directory}/host/{sample}.nohost.bam
	""".format(sample=sample, directory=directory)
	os.system(cmd)

# 比对
def mapping(sample, directory, rmhost):
	if rmhost:
		cmd = """
			mkdir -p {directory}/bam
			bwa mem -t 8 /home2/zhaowen/EVTpro/reference/EVTVP.fa \\
				{directory}/host/{sample}_1.fq.gz \\
				{directory}/host/{sample}_2.fq.gz \\
				| samtools view -bSh - | samtools sort -@ 8 - -o {directory}/bam/{sample}.bam
			samtools view -h -F 2052 {directory}/bam/{sample}.bam \\
				| samtools view -bSh - \\
				> {directory}/bam/{sample}.final.bam
			samtools index {directory}/bam/{sample}.final.bam
		""".format(sample=sample, directory=directory)
	else:
		cmd = """
			mkdir -p {directory}/bam
			bwa mem -t 8 /home2/zhaowen/EVTpro/reference/EVTVP.fa \\
				{directory}/cleandata/{sample}_1.fq.gz \\
				{directory}/cleandata/{sample}_2.fq.gz \\
				| samtools view -bSh - | samtools sort -@ 8 - -o {directory}/bam/{sample}.bam
			samtools view -h -F 2052 {directory}/bam/{sample}.bam \\
				| samtools view -bSh - \\
				> {directory}/bam/{sample}.final.bam
			samtools index {directory}/bam/{sample}.final.bam
		""".format(sample=sample, directory=directory)
	os.system(cmd)

# 直接统计比对结果
def typingStat(sample, directory, baseCover):
	bamFile = bamnostic.AlignmentFile(directory + "/bam/" + sample + ".final.bam", "rb")
	counter = []
	for read in bamFile:
		mapName = str(read).split("\t")[2].split("|")[1].split("_")[0]
		mapRegion = str(read).split("\t")[2].split("|")[1].split("_")[1]
		mapBaseCount = read.query_alignment_length
		if mapBaseCount <= baseCover:
			continue
		else:
			counter.append(mapName)

	outputDict = Counter(counter)
	output = sorted(outputDict.items(), key=lambda d: d[1])
	return output

# 使用在线blast校对结果
def ncbiBlastCheck(directory, sample):
	cmd = """
		mkdir -p {directory}/blast
		samtools fasta {directory}/bam/{sample}.final.bam > {directory}/blast/{sample}.fasta
	""".format(directory=directory, sample=sample)
	os.system(cmd)
	fastaFile = open(directory + "/blast/" + sample + ".fasta", "r")
	l = []
	for line in fastaFile:
		if line.startswith(">"):
			continue
		else:
			l.append(line)
	resultsDict = Counter(l)
	output = sorted(resultsDict.items(), key=lambda d: d[1])
	result_handle = NCBIWWW.qblast("blastn", "nt", output[0])
	ncbiOnlineCheckResult = open(directory + "/blast/" + sample + ".ncbiCheck.txt", "w")
	with open(directory + "/blast/" + sample + ".ncbiCheck.xml", "w") as save_to:
		save_to.write(result_handle.read())
		result_handle.close()
	tree = ET.ElementTree(file=directory + "/blast/" + sample + ".ncbiCheck.xml")
	hits = tree.iter("Hit")

	resultsDict = {}
	for hit in hits:
		for res in hit:
			if res.tag == "Hit_num":
				Hit_num = res.text
				resultsList = []

			if res.tag == "Hit_def":
				Hit_def = res.text
				resultsList.append(Hit_def)
			if res.tag == "Hit_accession":
				Hit_accession = res.text
				resultsList.append(Hit_accession)

			if res.tag == "Hit_hsps":
				for hsps in res:
					for hsp in hsps:
						if hsp.tag == "Hsp_bit-score":
							Hsp_bit_score = hsp.text
							resultsList.append(Hsp_bit_score)
						if hsp.tag == "Hsp_query-to":
							Hsp_query_to = hsp.text
							resultsList.append(Hsp_query_to)
						if hsp.tag == "Hsp_align-len":
							Hsp_align_len = hsp.text
							resultsList.append(Hsp_align_len)

							ncbiOnlineCheckResult.write("\t".join(resultsList) + "\n")

	ncbiOnlineCheckResult.close()

def main(sample, rawdata, directory, baseCover, rmhost, ncbicheck):
	filter(sample, rawdata, directory)
	if rmhost:
		mapToHost(sample, directory)
	mapping(sample, directory, rmhost)
	results = typingStat(sample, directory, baseCover)

	cmd = """
		mkdir {directory}/results
	""".format(directory=directory)
	os.system(cmd)

	resultsFile = open(directory + "/results/" + sample + ".type.txt", "w")

	print "#type\tcounts"
	resultsFile.write("#type\tcounts\n")
	for i in results:
		print i[0] + "\t" + str(i[1])
		resultsFile.write(i[0] + "\t" + str(i[1]) + "\n")

	if ncbicheck:
		ncbiBlastCheck(directory, sample)



if __name__ == "__main__":
	parser = argparse.ArgumentParser(description="EVT typing",
		prog="EVT.py",
		usage="python EVT.py -i <sampleID> -d <RawData Directory> -o <output Directory> -c <cut off>")
	parser.add_argument("-v", "--version", action="version",
		version="Version 2.0 20191120")
	parser.add_argument("-i", "--input", type=str,
		help="Input the sample ID")
	parser.add_argument("-d", "--dir", type=str,
		help="rawdata directory")
	parser.add_argument("-o", "--outputdir", type=str,
		help="output directory")
	parser.add_argument("-c", "--cutoff", type=int,
		help="baseCover cut off, default=70", default=70)
	parser.add_argument("-rmhost", "--removeHost", type=bool,
		help="remove hg19 reads, option", default=False)
	parser.add_argument("-ncbi", "--ncbiCheck", type=bool,
		help="ncbi blast check online, option", default=False)
	if len(sys.argv[1:]) == 0:
		parser.print_help()
		parser.exit()
	args = parser.parse_args()
	main(sample=args.input, rawdata=args.dir, directory=args.outputdir, baseCover=args.cutoff, rmhost=args.removeHost, ncbicheck=args.ncbiCheck)

