# coding=utf-8
# pzw
# 20190430
# v1.0

import sys
import os
import argparse
import pandas

# 参考组合并
def getReference(path):
	files = os.listdir(path)
	data = {}
	for file in files:
		name = file.split(".")[0]
		filex = open(path + "/" + file, "r")
		subID = {}
		n = 0
		for line in filex:
			if line.startswith("#"):
				continue
			else:
				counts = line.split("\t")[3]
				n += 1
				uniqueID = "ID_" + str(n).zfill(3)
				subID[uniqueID] = counts
		data[name] = subID
		filex.close()
	reference = pandas.DataFrame(data)
	return reference

def main(referencePath, outputFile):
	df = getReference(referencePath)
	df.to_csv(outputFile, sep="\t")

if __name__ == "__main__":
	parser = argparse.ArgumentParser(description="get SV reference database",
		prog="SV_reference_v1.0.py",
		usage="python SV_reference_v1.0.py -i <reference directory path> -o <output file>")
	parser.add_argument("-v", "--version", action="version", version="Version 1.0 20190430")
	parser.add_argument("-i", "--input", type=str, help="Input the path where the reference files locate in e.g:/path/reference/")
	parser.add_argument("-o", "--output", type=str, help="output the reference file e.g: reference.txt")
	if len(sys.argv[1:]) == 0:
		parser.print_help()
		parser.exit()
	args = parser.parse_args()
	main(referencePath=args.input, outputFile=args.output)