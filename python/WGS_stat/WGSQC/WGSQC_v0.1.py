# coding=utf-8
# pzw
# 20190514

import sys
import os
import argparse
import subprocess
import argcomplete

def main(function, option):

	now = os.path.dirname(sys.executable)

	function_dict = {
		"bedtools": "/home/zhaowen/workspace/software/WGSQC/function/bedtools2qc.py",
		"fastp": "/home/zhaowen/workspace/software/WGSQC/function/fastp2qc.py",
	}

	print("python " + os.path.join(now, function_dict[function]) + " " + " ".join(option))
	subprocess.call("python " + os.path.join(now, function_dict[function]) + " " + " ".join(option), shell=True)

if __name__ == "__main__":
	now = os.path.dirname(sys.executable)
	parser = argparse.ArgumentParser(description="WGS QC stats",
		prog="WGSQC_v0.1.py",
		usage="python WGSQC_v0.1.py [-h] <function> <function option>",
		formatter_class=argparse.RawTextHelpFormatter
	)
	parser.add_argument("-v", "--version", action="version", version="Version 0.1 20190514")
	parser.add_argument("function", choices=("bedtools", "fastp"),
		help='''
			bedtools                    统计bedtools结果文件
			fastp                       统计fastp结果文件
	''')
	parser.add_argument("option", nargs=argparse.REMAINDER, metavar="function option")
	argcomplete.autocomplete(parser)
	args = parser.parse_args()
	main(function=args.function, option=args.option)
