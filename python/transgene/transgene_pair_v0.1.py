# encoding:utf-8
# pzw
# 20190218

# v0.1
######################
import pandas as pd

#####################
countFilter = 10
gap = 100
result = pd.read_csv("result_test.txt", header=0, sep="\t")
outputFile = open("result_test_output.txt", "w")
#####################


## 过滤掉counts数达不到设定值的结果
resultAfterFilterCounts = result[result["counts"] > countFilter]

## 根据hostChr和hostPos列排序
resultAfterSort = resultAfterFilterCounts.sort_values(by=["hostChr", "hostPos", "plasmidPos"])
del resultAfterFilterCounts

## 重新定义index
resultAfterSort.index = range(len(resultAfterSort))

## 输出相邻的插入位点，获得pair
outputFile.write("hostChr\thostPos1\thostPos2\tplasmidPos1\tplasmidPos2\tfusionJunctionSeq1\tfusionJunctionSeq2\tcounts\n")
i = 0
while i < (len(resultAfterSort)-1):
	if resultAfterSort.loc[i+1]["hostChr"] == resultAfterSort.loc[i]["hostChr"]:
		if int(resultAfterSort.loc[i]["hostPos"]) > int(resultAfterSort.loc[i+1]["hostPos"]) - gap:
			if resultAfterSort.loc[i]["integrationType"] == resultAfterSort.loc[i+1]["integrationType"]:
				resultAfterSort.loc[i+1]["counts"] = resultAfterSort.loc[i]["counts"] + resultAfterSort.loc[i+1]["counts"]
				continue
			else:
				hostChr = resultAfterSort.loc[i]["hostChr"]
				hostPos1 = resultAfterSort.loc[i]["hostPos"]
				hostPos2 = resultAfterSort.loc[i+1]["hostPos"]
				plasmidPos1 = resultAfterSort.loc[i]["plasmidPos"]
				plasmidPos2 = resultAfterSort.loc[i+1]["plasmidPos"]
				fusionJunctionSeq1 = resultAfterSort.loc[i]["fusionJunctionSeq"]
				fusionJunctionSeq2 = resultAfterSort.loc[i+1]["fusionJunctionSeq"]
				counts = resultAfterSort.loc[i]["counts"] + resultAfterSort.loc[i+1]["counts"]
				outputList = [hostChr, str(hostPos1), str(hostPos2), str(plasmidPos1), str(plasmidPos2), fusionJunctionSeq1, fusionJunctionSeq2, str(counts)]
				outputFile.write("\t".join(outputList) + "\n")
	i += 1

outputFile.close()
del resultAfterSort
exit()
