# encoding:utf-8
# pzw
# 20190219
# python2.7
# v1.1 传参模式
# v1.0 替换为本地数据库
# v0.1
#######################
import pandas as pd
import sys

inputFile = open(sys.argv[1], "r")
outputFile = open(sys.argv[2], "w")
crigriDB = pd.read_csv(sys.argv[3], header=0, sep="\t")

for line in inputFile:
	if line.startswith("plasmidName"):
		print "print header"
		outputFile.write(line.split("\n")[0] + "\tGene stable ID\tTranscript stable ID\tChromosome/scaffold name\tGene start\tGene end\tGene description\tGene name\tGene Type\n")
	else:
		print "searching line"
		lineAfterSplit = line.split("\t")
		hostChr = lineAfterSplit[2]
		hostPos = lineAfterSplit[3]
		
		searchResult = crigriDB[(crigriDB["Chromosome/scaffold name"] == hostChr) & (crigriDB["Gene start"] <= int(hostPos)) & (crigriDB["Gene end"] >= int(hostPos))]
		result = searchResult.to_string(header=False, index=False, na_rep="") + "\n"
		if result.startswith("Empty"):
			print "can't find any match"
			outputFile.write(line)
		else:
			print result
			outputList = [line.split("\n")[0], "\t".join(result.split("  "))]
			outputFile.write("\t".join(outputList))

outputFile.close()
inputFile.close()
del crigriDB
print "task done"
exit()