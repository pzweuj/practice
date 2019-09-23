# coding=utf-8
# EVT pipeline
# pzw
# 20190923

import os
import sys
import bamnostic
from collections import Counter
import argparse

def filter(sample, rawdata, directory):
	cmd = """
		mkdir -p {directory}/QC
		mkdir -p {directory}/cleandata

		fastp -i {rawdata}/{sample}_R1_001.fastq.gz -I {rawdata}/{sample}_R2_001.fastq.gz \\
			-o {directory}/cleandata/{sample}_1.fq.gz -O {directory}/cleandata/{sample}_2.fq.gz \\
			-w 8 -j {directory}/QC/{sample}.json -h {directory}/QC/{sample}.html
	""".format(sample=sample, rawdata=rawdata, directory=directory)
	os.system(cmd)

def mapping(sample, directory):
	cmd = """
		mkdir -p {directory}/bam
		bwa mem -t 8 /home2/zhaowen/project/EVTpro/reference/EVTREF_VP.fa \\
			{directory}/cleandata/{sample}_1.fq.gz \\
			{directory}/cleandata/{sample}_2.fq.gz \\
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
	mapping(sample, directory)
	results = typingStat(sample, directory, baseCover)
	print "#type\tcounts"
	for i in results:
		print i[0] + "\t" + str(i[1])


if __name__ == "__main__":
	parser = argparse.ArgumentParser(description="EVT typing",
		prog="EVT.py",
		usage="python EVT.py -i <sampleID> -d <RawData Directory> -o <output Directory> -c <cut off>")
	parser.add_argument("-v", "--version", action="version",
		version="Version 0.1 20190923")
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

