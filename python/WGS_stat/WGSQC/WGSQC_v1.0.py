# coding=utf-8
# pzw
# 20190531

import sys
import os
import argparse
import subprocess
import argcomplete

def main(function, option):

	now = os.path.abspath(os.path.dirname(sys.argv[0]))

	function_dict = {
		"bedtools": "/function/bedtools2qc.py",
		"fastp": "/function/fastp2qc.py",
	}

	print("python " + now + function_dict[function] + " " + " ".join(option))
	subprocess.call("python " + now + function_dict[function] + " " + " ".join(option), shell=True)

if __name__ == "__main__":
	now = os.path.dirname(sys.executable)
	parser = argparse.ArgumentParser(
		description="WGS QC stats  @PZW",
		prog="WGSQC_v1.0.py",
		usage="python WGSQC_v1.0.py [-h] <function> <function option>",
		formatter_class=argparse.RawTextHelpFormatter
	)
	parser.add_argument("-v", "--version", action="version", version="Version 1.0 20190531")
	parser.add_argument("function", choices=("bedtools", "fastp"),
		help='''
			bedtools                    analysis bedtools genenocov results
			fastp                       analysis fastp json file
	''')
	parser.add_argument("option", nargs=argparse.REMAINDER, metavar="function option")
	argcomplete.autocomplete(parser)
	if len(sys.argv[1:]) == 0:
		parser.print_help()
		parser.exit()
	args = parser.parse_args()
	main(function=args.function, option=args.option)
