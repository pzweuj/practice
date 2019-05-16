# pzw
# 20190516

import sys
import argparse

def main(resultsFile, bed, outputFile):
	results = open(resultsFile, "r")
	exonbed = open(bed, "r")
	output = open(outputFile, "w")

	resultsDict = {}
	i = 1
	for line in exonbed:

		if line.startswith("#"):
			continue
		else:
			lineAfterSplit = line.split("\t")
			chrom = lineAfterSplit[0]
			start = int(lineAfterSplit[1])
			end = int(lineAfterSplit[2])
			gene = lineAfterSplit[3].split("\n")[0]

		resultsDict[str(i)] = [chrom, start, end, gene]
		i += 1
	exonbed.close()

	for line2 in results:
		l = []
		if line2.startswith("#"):
			output.write(line2.split("\n")[0] + "\texonLengthPercent\n")
		else:
			lineAfterSplit2 = line2.split("\t")
			chrom2 = lineAfterSplit2[0]
			start2 = int(lineAfterSplit2[1])
			end2 = int(lineAfterSplit2[2])
			gene2 = lineAfterSplit2[3].split("\n")[0]
			length_all = end2 - start2

			exonLength = 0
			for m in resultsDict:
				length_tmp = 0
				if resultsDict[m][0] == chrom2:
					if resultsDict[m][1] < start2 and start2 < resultsDict[m][2]:
						length_tmp = resultsDict[m][2] - start2
					elif start2 <= resultsDict[m][1] and resultsDict[m][2] <= end2:
						length_tmp = resultsDict[m][2] - resultsDict[m][1]
					elif resultsDict[m][1] <= start2 and end2 <= resultsDict[m][2]:
						length_tmp = end2 - start2
					elif resultsDict[m][1] < end2 and end2 < resultsDict[m][2]:
						length_tmp = end2 - resultsDict[m][1]
					else:
						continue
				exonLength += length_tmp
			print(line2, exonLength)

			exonprecent = float(exonLength) / float(length_all)

			output.write(line2.split("\n")[0] + "\t" + str(exonprecent) + "\n")

	output.close()
	results.close()

if __name__ == "__main__":
	parser = argparse.ArgumentParser(description="get exon content",
		prog="BRCA_add_exon.py",
		usage="python3 BRCA_add_exon.py -i <input file> -o <output file> -b <bedfile>")
	parser.add_argument("-v", "--version", action="version", version="Version 0.1 20190516")
	parser.add_argument("-i", "--input", type=str, help="Input the file from BRCA summary")
	parser.add_argument("-o", "--output", type=str, help="Output the result file")
	parser.add_argument("-b", "--bedfile", type=str, default="/home/zhaowen/workspace/software/BRCATK/function/hg19.exons.BRCA.bed", help="Input the bedfile which was used to get exon. default: hg19.exons.BRCA.bed")
	if len(sys.argv[1:]) == 0:
		parser.print_help()
		parser.exit()
	args = parser.parse_args()
	main(resultsFile=args.input, bed=args.bedfile, outputFile=args.output)


