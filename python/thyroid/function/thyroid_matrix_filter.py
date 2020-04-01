# coding=utf-8
# pzw
# 20200325

import sys
import argparse

def matrix_filter(matrixFile, cleanMatrix):
	empty = []
	thyroid_matrix = open(matrixFile, "r")
	data_dict = {}
	important_gene = ["AKT1", "ALK", "BRAF", "CTNNB1", "GNAS", "HRAS", "KRAS", "PIK3CA", "PTEN", "RET", "TERT", "TP53", "TSHR"]
	for i in thyroid_matrix:
		if i.startswith("Gene"):
			GeneList = i.split("\t")
			for g in range(len(GeneList)):
				if g != 0:
					if GeneList[g] not in important_gene:
						empty.append(g - 1)
		elif not (i.startswith("H") or i.startswith("T") or i.startswith("U")):
			continue
		else:
			xx = i.split("\n")[0].split("\t")
			variantCount = len(xx) - 1
			ID = xx[0]
			data_dict[ID] = xx[1:]
	thyroid_matrix.close()

	j = 0
	while j < variantCount:
		seq = ""
		for k in data_dict.keys():
			seq = seq + data_dict[k][j]
		if seq == "":
			empty.append(j)
		j += 1

	empty = sorted(list(set(empty)), key=int)


	thyroid_matrix = open(matrixFile, "r")
	thyroid_matrix_filter = open(cleanMatrix, "w")
	for line in thyroid_matrix:
		ls = line.split("\n")[0].split("\t")
		output = []
		for i in range(len(ls)):
			if (i - 1) not in empty:
				output.append(ls[i])
		thyroid_matrix_filter.write("\t".join(output) + "\n")

	thyroid_matrix.close()
	thyroid_matrix_filter.close()


def main(matrixFile, cleanMatrix):
	matrix_filter(matrixFile, cleanMatrix)

if __name__ == "__main__":
	parser = argparse.ArgumentParser(description="Matrix Filter",
		prog="thyroid_matrix_filter.py",
		usage="python3 thyroid_matrix_filter.py -i <raw matrix file> -o <clean matrix file>")
	group = parser.add_mutually_exclusive_group()
	parser.add_argument("-v", "--version", action="version",
		version="Version 0.1 20200330")
	parser.add_argument("-i", "--input", type=str,
		help="Input the matrix file")
	parser.add_argument("-o", "--output", type=str,
		help="output the clean matrix")
	if len(sys.argv[1:]) == 0:
		parser.print_help()
		parser.exit()
	args = parser.parse_args()
	main(matrixFile=args.input, cleanMatrix=args.output)
