# coding=utf-8
# pzw
# 20190517

import sys
import os
import argparse

def getChromCoverage(inputfile, chrom="all"):
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

		if chrom == "all":
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


def readChromName(inputfile):
	inputFile = open(inputfile, "r")
	l = []
	for line in inputFile:
		lineAS = line.split("\t")
		chromosome = lineAS[0]
		l.append(chromosome)
	news_l = list(set(l))
	news_l.sort(key=l.index)
	inputFile.close()


	finalList = []
	for i in range(len(news_l)):
		if len(news_l[i].split("_")) <= 1:
			finalList.append(news_l[i])

	return finalList



def main(inputFile, outputChr, outputFile):

	if outputFile:

		output = open(outputFile, "w")

		chromList = readChromName(inputFile)

		coverage_list = ["Coverage"]
		depth_cov_list = ["CoverRegionDepth"]
		depth_all_list = ["Depth"]

		for m in chromList:
			chromi = m.split("\n")[0]
			chromeResults = getChromCoverage(inputFile, chromi)
			print chromi, chromeResults
			coverage_list.append(chromeResults[0])
			depth_cov_list.append(chromeResults[1])
			depth_all_list.append(chromeResults[2])

		summaryResults = getChromCoverage(inputFile, "all")
		print "summary", summaryResults
		chromList.append("summary")
		coverage_list.append(summaryResults[0])
		depth_cov_list.append(summaryResults[1])
		depth_all_list.append(summaryResults[2])

		output.write("\t" + "\t".join(chromList) + "\n")
		output.write("\t".join(coverage_list) + "\n")
		output.write("\t".join(depth_cov_list) + "\n")
		output.write("\t".join(depth_all_list) + "\n")

		output.close()

	if outputChr:
		print "chromosome\t[coverage, mapped depth, depth]"
		chromList = outputChr.split(",")
		for m in chromList:
			chromi = m.split("\n")[0]
			print chromi, "\t", getChromCoverage(inputFile, chromi)


if __name__ == "__main__":
	parser = argparse.ArgumentParser(description="analysis bedtools output",
		prog="bedtools2qc.py",
		usage="python bedtools2qc.py -i <bedtools.output> -o <results>")
	group = parser.add_mutually_exclusive_group()
	parser.add_argument("-v", "--version", action="version",
		version="Version 1.5 20190520")
	parser.add_argument("-i", "--input", type=str,
		help="Input the file which output from 'bedtools genomecov -ibam bam -g reference -bga'")
	group.add_argument("-o", "--output", type=str,
		help="output the results")
	group.add_argument("-chr", "--chromosome", type=str,
		help="output select chromosome stats e.g1: -chr chr1, e.g2: -chr chr1,chr3,chr5, e.g3: -chr all")
	if len(sys.argv[1:]) == 0:
		parser.print_help()
		parser.exit()
	args = parser.parse_args()
	if args.output and args.chromosome:
		sys.exit("-o and -chr must be only one!")
	main(inputFile=args.input, outputChr=args.chromosome, outputFile=args.output)
