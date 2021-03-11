# coding=utf-8
# pzw
# 20210311

import os
import sys

a = open("a.txt", "r")
b = open("results.txt", "w")
for line in a:
	l = line.replace("\n", "").split("\t")
	chrom = "chr" + l[0]
	start = l[1]
	end = l[2]

	cmd = """
		samtools faidx /home/bioinfo/ubuntu/database/hg19/ucsc.hg19.fasta \\
			{chrom}:{start}-{end}
	""".format(chrom=chrom, start=start, end=end)
	results = os.popen(cmd)
	res = results.read()
	check = res.split("\n")[0]
	output = "".join(res.split("\n")[1:]).upper()

	b.write("\t".join([check, output]) + "\n")
