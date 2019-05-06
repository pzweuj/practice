# coding=utf-8
# pzw
# 20190506
# v0.1

import subprocess
import argparse
import sys

def main(inputfile, sampleID, method, outputDir):
	if method == "DQ":
		subprocess.call(["Rscript", "function/SV_plot_DQ.R", inputfile, sampleID, outputDir])
	elif method == "Z":
		subprocess.call(["Rscript", "function/SV_plot_Z.R", inputfile, sampleID, outputDir])
	else:
		print("please enter method (DQ or Z)")

if __name__ == "__main__":
	parser = argparse.ArgumentParser(description="plot function",
		prog="SV_plot_v0.1.py",
		usage="python3 SV_plot_v0.1.py -i <input file> -I <ID> -m <method> -o <output directory>")
	parser.add_argument("-v", "--version", action="version", version="Version 0.1 20190506")
	parser.add_argument("-i", "--input", type=str, help="Input the the file from SV_analysis")
	parser.add_argument("-I", "--ID", type=str, help="Input the test sample ID e.g: XXX")
	parser.add_argument("-m", "--method", type=str, help="chose a method e.g: Z or DQ")
	parser.add_argument("-o", "--output", type=str, help="output directory e.g: /path/to/results")
	if len(sys.argv[1:]) == 0:
		parser.print_help()
		parser.exit()
	args = parser.parse_args()
	main(inputfile=args.input, sampleID=args.ID, method=args.method, outputDir=args.output)
