# pzw
# 20190514

import sys
import json
import os
import argparse

def main(jsonFile, outputFile):
	QC_json = open(jsonFile, "r")
	output = open(outputFile, "w")
	jsonFile = json.load(QC_json)

	beforeFilter = jsonFile["summary"]["before_filtering"]
	afterFilter = jsonFile["summary"]["after_filtering"]

	before_reads = beforeFilter["total_reads"]
	before_bases = beforeFilter["total_bases"]
	before_GC = "%.2f" % (beforeFilter["gc_content"] * 100) + "%"
	before_Q20 = "%.2f" % (beforeFilter["q20_rate"] * 100) + "%"
	before_Q30 = "%.2f" % (beforeFilter["q30_rate"] * 100) + "%"
	duplicateRate = "%.2f" % (jsonFile["duplication"]["rate"] * 100) + "%"

	after_reads = afterFilter["total_reads"]
	after_bases = afterFilter["total_bases"]
	after_GC = "%.2f" % (afterFilter["gc_content"] * 100) + "%"
	after_Q20 = "%.2f" % (afterFilter["q20_rate"] * 100) + "%"
	after_Q30 = "%.2f" % (afterFilter["q30_rate"] * 100) + "%"

	print "raw reads: ", before_reads
	print "raw bases: ", before_bases
	print "raw GC content: ", before_GC
	print "raw Q20: ", before_Q20
	print "raw Q30: ", before_Q30
	print "duplication: ", duplicateRate
	print "------------------------------------------"
	print "clean reads: ", after_reads
	print "clean bases: ", after_bases
	print "clean GC content: ", after_GC
	print "clean Q20: ", after_Q20
	print "clean Q30: ", after_Q30

	output.write("raw reads: " + str(before_reads) + "\n")
	output.write("raw bases: " + str(before_bases) + "\n")
	output.write("raw GC content: " + str(before_GC) + "\n")
	output.write("raw Q20: " + before_Q20 + "\n")
	output.write("raw Q30: " + before_Q30 + "\n")
	output.write("duplication: " + duplicateRate + "\n")
	output.write("------------------------------------------\n")
	output.write("clean reads: " + str(after_reads) + "\n")
	output.write("clean bases: " + str(after_bases) + "\n")
	output.write("clean GC content: " + after_GC + "\n")
	output.write("clean Q20: " + after_Q20 + "\n")
	output.write("clean Q30: " + after_Q30 + "\n")

	output.close()
	QC_json.close()

	print "fastp output file analysis done"

if __name__ == "__main__":
	parser = argparse.ArgumentParser(description="analysis fastp output",
		prog="fastp2qc.py",
		usage="python fastp2qc.py -i <fastp.json> -o <results>")
	parser.add_argument("-v", "--version", action="version",
		version="Version 1.0 20190514")
	parser.add_argument("-i", "--input", type=str,
		help="Input the file which output 'fastp.json'")
	parser.add_argument("-o", "--output", type=str,
		help="output the results")
	if len(sys.argv[1:]) == 0:
		parser.print_help()
		parser.exit()
	args = parser.parse_args()
	main(jsonFile=args.input, outputFile=args.output)