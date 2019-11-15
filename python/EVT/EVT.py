# coding=utf-8
# EVT pipeline
# pzw
# 20191114

import os
import sys
import bamnostic
from collections import Counter
import argparse
import json
import commands

def filter(sample, rawdata, directory):
	cmd = """
		mkdir -p {directory}/QC
		mkdir -p {directory}/cleandata

		fastp -i {rawdata}/{sample}_R1_001.fastq.gz -I {rawdata}/{sample}_R2_001.fastq.gz \\
			-o {directory}/cleandata/{sample}_1.fq.gz -O {directory}/cleandata/{sample}_2.fq.gz \\
			-w 8 -j {directory}/QC/{sample}.json -h {directory}/QC/{sample}.html
	""".format(sample=sample, rawdata=rawdata, directory=directory)
	os.system(cmd)

def mapToHost(sample, directory):
	cmd = """
		mkdir -p {directory}/filter
		bwa mem -t 8 /home/zhaowen/workspace/database/human/hg19.fa \\
			{directory}/cleandata/{sample}_1.fq.gz \\
			{directory}/cleandata/{sample}_2.fq.gz \\
			> {directory}/filter/{sample}.host.sam
		samtools view {directory}/filter/{sample}.host.sam -F 2052 | wc -l > {directory}/filter/{sample}.hostCount.txt
		samtools view -bSh {directory}/filter/{sample}.host.sam -f 12 -F 256 > {directory}/filter/{sample}.nohost.bam
		rm {directory}/filter/{sample}.host.sam
		samtools sort -n {directory}/filter/{sample}.nohost.bam -@ 8 -o {directory}/filter/{sample}.nohost.sorted.bam
		bedtools bamtofastq -i {directory}/filter/{sample}.nohost.sorted.bam \\
			-fq {directory}/filter/{sample}.rmhost_1.fastq \\
			-fq2 {directory}/filter/{sample}.rmhost_2.fastq
		gzip {directory}/filter/{sample}.rmhost_1.fastq
		gzip {directory}/filter/{sample}.rmhost_2.fastq
	""".format(sample=sample, directory=directory)
	os.system(cmd)


def mapping(sample, directory):
	cmd = """
		mkdir -p {directory}/bam
		bwa mem -t 8 /home/zhaowen/workspace/database/EVT/EVTREF_VP.fa \\
			{directory}/filter/{sample}.rmhost_1.fastq.gz \\
			{directory}/filter/{sample}.rmhost_2.fastq.gz \\
			| samtools view -bSh - | samtools sort -@ 8 - -o {directory}/bam/{sample}.bam
		samtools view {directory}/bam/{sample}.bam -H > {directory}/bam/{sample}.header
		samtools view {directory}/bam/{sample}.bam -F 2052 -bSh > {directory}/bam/{sample}.temp.bam
		samtools index {directory}/bam/{sample}.temp.bam
	""".format(sample=sample, directory=directory)
	os.system(cmd)
	bamFile = bamnostic.AlignmentFile(directory + "/bam/" + sample + ".temp.bam", "rb")
	finalSam = open(directory + "/bam/" + sample + ".sam", "w")
	for read in bamFile:
		if str(read).split("\t")[6] == "=":
			finalSam.write(str(read) + "\n")
	finalSam.close()

	cmd2 = """
		cat {directory}/bam/{sample}.header {directory}/bam/{sample}.sam | samtools view -bSh - > {directory}/bam/{sample}.final.bam
		rm {directory}/bam/{sample}.header {directory}/bam/{sample}.sam {directory}/bam/{sample}.temp.bam*
		samtools index {directory}/bam/{sample}.final.bam
	""".format(sample=sample, directory=directory)
	os.system(cmd2)

def typingStat(sample, directory, baseCover):
	bamFile = bamnostic.AlignmentFile(directory + "/bam/" + sample + ".final.bam", "rb")
	counter = []
	for read in bamFile:
		mapName = str(read).split("\t")[2].split("|")[1]
		mapBaseCount = read.query_alignment_length
		if mapBaseCount <= baseCover:
			continue
		else:
			counter.append(mapName)

	outputDict = Counter(counter)
	output = sorted(outputDict.items(), key=lambda d: d[1])
	return output


def main(sample, rawdata, directory, baseCover):
	filter(sample, rawdata, directory)
	mapToHost(sample, directory)
	mapping(sample, directory)
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


if __name__ == "__main__":
	parser = argparse.ArgumentParser(description="EVT typing",
		prog="EVT.py",
		usage="python EVT.py -i <sampleID> -d <RawData Directory> -o <output Directory> -c <cut off>")
	parser.add_argument("-v", "--version", action="version",
		version="Version 0.2 20191114")
	parser.add_argument("-i", "--input", type=str,
		help="Input the sample ID")
	parser.add_argument("-d", "--dir", type=str,
		help="rawdata directory")
	parser.add_argument("-o", "--outputdir", type=str,
		help="output directory")
	parser.add_argument("-c", "--cutoff", type=int,
		help="baseCover cut off", default=70)
	if len(sys.argv[1:]) == 0:
		parser.print_help()
		parser.exit()
	args = parser.parse_args()
	main(sample=args.input, rawdata=args.dir, directory=args.outputdir, baseCover=args.cutoff)

