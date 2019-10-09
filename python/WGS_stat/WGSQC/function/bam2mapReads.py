# pzw
# 20191009

import commands
import os
import sys
import argparse

def main(bamFile):
	cmd = "samtools view -F 2052 " + bamFile + " | wc -l"
	cmdOutput = commands.getstatusoutput(cmd)
	mapReads = cmdOutput[1]
	print mapReads

if __name__ == "__main__":
	parser = argparse.ArgumentParser(description="analysis bamFile mapping reads",
		prog="bam2mapReads.py",
		usage="python bam2mapReads.py -i <bamFile>")
	parser.add_argument("-v", "--version", action="version",
		version="Version 0.1 20191009")
	parser.add_argument("-i", "--input", type=str,
		help="Input bam file")
	if len(sys.argv[1:]) == 0:
		parser.print_help()
		parser.exit()
	args = parser.parse_args()
	main(bamFile=args.input)