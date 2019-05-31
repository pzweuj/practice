# coding=utf-8
# pzw
# 20190529
# 过滤sam文件中比对到多个位置以及错配大于1个碱基的reads

from cigar import Cigar

samFile = open("test.sam", "r")
filterDict = {}
for line in samFile:
	if line.startswith("@"):
		continue
	else:
		l = line.split("\t")
		readID = l[0]
		chrName = l[2]
		filterDict[readID] = []
		filterDict[readID].append(chrName)
		filterDict[readID] = set(filterDict[readID])

samFile.close()

samFile2 = open("QH046cDNA.sam", "r")
for line2 in samFile2:
	if line2.startswith("@"):
		print line2
	else:
		l2 = line2.split("\t")
		readID2 = l[0]
		if len(filterDict[readID2]) > 2:
			continue
		else:
			seqLength = len(l2[9])
			cigar = Cigar(l2[5])
			cigarList = list(cigar.items())
			mapp = 0
			for i in cigarList:
				if i[1] == "M":
					mapp += i[0]
				else:
					continue
			if seqLength - mapp <= 1:
				print line2
samFile2.close()