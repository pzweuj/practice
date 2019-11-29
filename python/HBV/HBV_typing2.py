# coding=utf-8
# pzw
# 20191129
# v2.0

import sys
import re
import argparse

def CreateReadDict(blastResultFile):
	blastResults = open(blastResultFile, "r")
	i = 0
	n = 0
	readNameList = []
	bestHit = []
	readCountsDict = {}

	for line in blastResults:

		# 跳过上面14行
		if i <= 14:
			i += 1
			continue
		else:
			if line.startswith("Query= "):
				readName = line.split("Query= ")[1]
				readNameList.append(readName)
				readCounts = int(line.split("_counts_")[1])
				readCountsDict[readName] = readCounts
				n = 0
			else:
				n += 1

				# 提取Query下第6行的best hit
				if n == 6:
					a = re.split(" ", line.split("\n")[0])
					while "" in a:
						a.remove("")
					bestHit.append(a)
	blastResults.close()
	hitResults = zip(readNameList, bestHit)
	hitReads = {"A": 0, "B": 0, "C": 0, "D": 0, "E": 0, "F": 0, "G": 0, "H": 0}

	# 统计
	nReads = 0
	for result in hitResults:
		nReads += readCountsDict[result[0]]

		# 避免出现No hit found的情况
		if len(result[1]) >= 3:
			genotype = result[1][0].split("|")[1]
			hitReads[genotype] += readCountsDict[result[0]]

	hitReadSort = sorted(hitReads.items(), lambda x, y: cmp(x[1], y[1]), reverse=True)

	return hitReadSort, nReads



def main(blastResultFile, readCut):
	ReadDict = CreateReadDict(blastResultFile)
	nReads = ReadDict[1]
	ReadDictResults = ReadDict[0]

	bestHitPercent = float(ReadDictResults[0][1]) / nReads
	if nReads >= readCut:
		if bestHitPercent >= 0.85:
			print ReadDictResults[0][0] + "\t" + str("%.2f" % float(bestHitPercent*100)) + "%\t" + str(ReadDictResults[0][1])
		elif bestHitPercent >= 0.45:
			print ReadDictResults[0][0] + "\t" + str("%.2f" % float(bestHitPercent*100)) + "%\t" + str(ReadDictResults[0][1])
			print ReadDictResults[1][0] + "\t" + str("%.2f" % (float(ReadDictResults[1][1]) / nReads * 100)) + "%\t" + str(ReadDictResults[1][1])
		elif bestHitPercent >= 0.3:
			print ReadDictResults[0][0] + "\t" + str("%.2f" % float(bestHitPercent*100)) + "%\t" + str(ReadDictResults[0][1])
			print ReadDictResults[1][0] + "\t" + str("%.2f" % (float(ReadDictResults[1][1]) / nReads * 100)) + "%\t" + str(ReadDictResults[1][1])
			print ReadDictResults[1][0] + "\t" + str("%.2f" % (float(ReadDictResults[2][1]) / nReads * 100)) + "%\t" + str(ReadDictResults[2][1])
		else:
			print "未找到任何分型\t-"
	else:
		print "未达分型阈值，未找到任何分型\t-"


if __name__ == "__main__":
	parser = argparse.ArgumentParser(description="HBV typing @PZW",
		prog="HBV_typing.py",
		usage="python HBV_typing.py -i <blast output> -c <cut off>")
	parser.add_argument("-v", "--version", action="version",
		version="Version 2.0 20191129")
	parser.add_argument("-i", "--input", type=str,
		help="Input the blast result file")
	parser.add_argument("-c", "--cutoff", type=int,
		help="min total reads amount")
	if len(sys.argv[1:]) == 0:
		parser.print_help()
		parser.exit()
	args = parser.parse_args()
	main(blastResultFile=args.input, readCut=args.cutoff)
