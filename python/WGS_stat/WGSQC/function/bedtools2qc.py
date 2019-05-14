# coding=utf-8
# pzw
# 20190513

import sys
import os
import argparse

def getChromCoverage(inputfile, chrom):
	inputFile = open(inputfile, "r")
	length = 0
	covLength = 0
	counts = 0
	for line in inputFile:
		lineAS = line.split("\t")
		chromosome = lineAS[0]
		start = lineAS[1]
		end = lineAS[2]
		count = lineAS[3].split("\n")[0]
		if chromosome == chrom:
			length_tmp = int(end) - int(start)
			length += length_tmp
			counts += (int(count) * length_tmp)
			if count != "0":
				covLength_tmp = length_tmp
				covLength += covLength_tmp

	coverage = "%.2f" % ((float(covLength) / float(length)) * 100) + "%"
	depth_cov = "%.2f" % (float(counts) / float(covLength))
	depth_all = "%.2f" % (float(counts) / float(length))
	inputFile.close()

	return [coverage, depth_cov, depth_all]


def main(inputFile, outputFile):

	output = open(outputFile, "w")

	chromList = ["chrX", "chrY", "chrM"]
	i = 1
	while i <= 22:
		chrom = "chr" + str(i)
		i += 1
		chromList.append(chrom)


	coverage_list = ["覆盖度"]
	depth_cov_list = ["覆盖区域的平均深度"]
	depth_all_list = ["WGS平均深度"]

	for m in chromList:
		chromi = m.split("\n")[0]
		print chromi, getChromCoverage(inputFile, chromi)
		coverage_list.append(getChromCoverage(inputFile, chromi)[0])
		depth_cov_list.append(getChromCoverage(inputFile, chromi)[1])
		depth_all_list.append(getChromCoverage(inputFile, chromi)[2])

	output.write("\t" + "\t".join(chromList) + "\n")
	output.write("\t".join(coverage_list) + "\n")
	output.write("\t".join(depth_cov_list) + "\n")
	output.write("\t".join(depth_all_list) + "\n")

	output.close()

if __name__ == "__main__":
	parser = argparse.ArgumentParser(description="analysis bedtools output",
		prog="bedtools2qc.py",
		usage="python bedtools2qc.py -i <bedtools.output> -o <results>")
	parser.add_argument("-v", "--version", action="version",
		version="Version 1.1 20190509")
	parser.add_argument("-i", "--input", type=str,
		help="Input the file which output from 'bedtools genomecov -ibam bam -g reference -bga'")
	parser.add_argument("-o", "--output", type=str,
		help="output the results")
	if len(sys.argv[1:]) == 0:
		parser.print_help()
		parser.exit()
	args = parser.parse_args()
	main(inputFile=args.input, outputFile=args.output)
