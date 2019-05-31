# pzw
# 20190523

import sys
import json
import os
import argparse

def main(jsonfile, outputFile):
	QC_json = open(jsonfile, "r")
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

	print("raw reads: ", before_reads)
	print("raw bases: ", before_bases)
	print("raw GC content: ", before_GC)
	print("raw Q20: ", before_Q20)
	print("raw Q30: ", before_Q30)
	print("duplication: ", duplicateRate)
	print("------------------------------------------")
	print("clean reads: ", after_reads)
	print("clean bases: ", after_bases)
	print("clean GC content: ", after_GC)
	print("clean Q20: ", after_Q20)
	print("clean Q30: ", after_Q30)

	jsonName = ""
	jsonAS = jsonfile.split("/")
	for jsons in jsonAS:
		if jsons.__contains__("json"):
			jsonName = jsons.split(".json")[0]
	output.write(jsonName + "\n")
	output.write("rawdata_reads\trawdata_bases\trawdata_dups\trawdata_GC\trawdata_Q20\trawdata_Q30\tcleandata_reads\tcleandata_bases\tcleandata_GC\tcleandata_Q20\tcleandata_Q30\n")

	outputList = [str(before_reads), str(before_bases), duplicateRate, str(before_GC), before_Q20, before_Q30, str(after_reads), str(after_bases), after_GC, after_Q20, after_Q30]
	output.write("\t".join(outputList) + "\n")

	output.close()
	QC_json.close()

	print("fastp output file analysis done")

if __name__ == "__main__":
	parser = argparse.ArgumentParser(description="analysis fastp output",
		prog="fastp2qc.py",
		usage="python fastp2qc.py -i <fastp.json> -o <results>")
	parser.add_argument("-v", "--version", action="version",
		version="Version 1.1 20190523")
	parser.add_argument("-i", "--input", type=str,
		help="Input the file which output 'fastp.json'")
	parser.add_argument("-o", "--output", type=str,
		help="output the results")
	if len(sys.argv[1:]) == 0:
		parser.print_help()
		parser.exit()
	args = parser.parse_args()
	main(jsonfile=args.input, outputFile=args.output)