# encoding:utf-8
# pzw
# 20190219
# python2.7

# update log
# v0.5 修改为传参模式
# v0.4 添加运行log
# v0.3 修复warning
# v0.2 修复一些bug
# v0.1
######################
import pandas as pd
import time
import sys
#####################
result = pd.read_csv(sys.argv[1], header=0, sep="\t")
outputFile = open(sys.argv[2], "w")
countFilter = int(sys.argv[3]) 
gap = int(sys.argv[4])
#####################

####### test ########
# result = pd.read_csv("result_test.txt", header=0, sep="\t")
# outputFile = open("result_test_output.txt", "w")
# countFilter = 10
# gap = 100
#####################


start = time.clock()
## 过滤掉counts数达不到设定值的结果
print "setting gap is", gap
print "extract counts >", countFilter
resultAfterFilterCounts = result[result["counts"] > countFilter]
plasmidName = result["plasmidName"][0]

## 根据hostChr和hostPos列排序
resultAfterSort = resultAfterFilterCounts.sort_values(by=["hostChr", "hostPos", "plasmidPos"])
print "sorting done"
resultAfterSort.to_csv("results.sort.txt", sep="\t", index=False)
print "output sorted file"
del resultAfterFilterCounts

## 重新定义index
resultAfterSort.index = range(len(resultAfterSort))

## 输出相邻的插入位点，获得pair
outputFile.write("hostChr\thostPos1\thostPos2\tplasmidPos1\tplasmidPos2\tfusionJunctionSeq1\tfusionJunctionSeq2\tcounts\n")

def getDirect(x, y):
	if resultAfterSort.loc[x]["integrationType"].split("/")[0] == y:
		return "right"
	else:
		return "left"

## 合并插入位点相近的结果读数
i = 0
while i < (len(resultAfterSort)-1):
	## 当方向为宿主/质粒
	if getDirect(i, plasmidName) == "left":
		## 插入方式一致证明是同向且同染色体
		if resultAfterSort.loc[i+1]["integrationType"] == resultAfterSort.loc[i]["integrationType"]:
			## 位置相隔小于阈值
			if resultAfterSort.loc[i]["hostPos"] > resultAfterSort.loc[i+1]["hostPos"] - gap:
				## 以pos比较大的为结果
				resultAfterSort.iat[i+1, 6] = resultAfterSort.loc[i]["counts"] + resultAfterSort.loc[i+1]["counts"]
	## 当方向为质粒/宿主
	if getDirect(i, plasmidName) == "right":
		## 插入方式一致证明是同向且同染色体
		if resultAfterSort.loc[i+1]["integrationType"] == resultAfterSort.loc[i]["integrationType"]:
			## 位置相隔小于阈值
			if resultAfterSort.loc[i]["hostPos"] > resultAfterSort.loc[i+1]["hostPos"] - gap:
				## 以pos比较小的为结果
				resultAfterSort.iat[i, 6] = resultAfterSort.loc[i]["counts"] + resultAfterSort.loc[i+1]["counts"]
	i += 1
print "combine results' counts"
	
## 输出结果
j = 0
while j < (len(resultAfterSort)-1):
	## 判断方向是否相反
	if getDirect(j, plasmidName) == "left" and getDirect(j+1, plasmidName) == "right":
		print j, j+1, "opposite direction, pass"
		## 判断是否在同一染色体上
		if resultAfterSort.loc[j+1]["hostChr"] == resultAfterSort.loc[j]["hostChr"]:
			print j, j+1, "same chrom/contig, pass"
			## 判断间隔是否在阈值之内
			if resultAfterSort.loc[j]["hostPos"] > resultAfterSort.loc[j+1]["hostPos"] - gap:
				print j, j+1, "gap <", gap, "pass, output to results"
				hostChr = resultAfterSort.loc[j]["hostChr"]
				hostPos1 = resultAfterSort.loc[j]["hostPos"]
				hostPos2 = resultAfterSort.loc[j+1]["hostPos"]
				plasmidPos1 = resultAfterSort.loc[j]["plasmidPos"]
				plasmidPos2 = resultAfterSort.loc[j+1]["plasmidPos"]
				fusionJunctionSeq1 = resultAfterSort.loc[j]["fusionJunctionSeq"]
				fusionJunctionSeq2 = resultAfterSort.loc[j+1]["fusionJunctionSeq"]
				counts = str(resultAfterSort.loc[j]["counts"] + resultAfterSort.loc[j+1]["counts"]) + "/" + str(resultAfterSort.loc[j]["counts"]) + "/" + str(resultAfterSort.loc[j+1]["counts"])
				outputList = [hostChr, str(hostPos1), str(hostPos2), str(plasmidPos1), str(plasmidPos2), fusionJunctionSeq1, fusionJunctionSeq2, counts]
				outputFile.write("\t".join(outputList) + "\n")
				j += 1
			else:
				print j, j+1, "gap >", gap, "fail"
				j += 1
				continue
		else:
			print j, j+1, "different chrom/contig, fail"
			j += 1
			continue
	else:
		print j, j+1, "same direction, fail"
		j += 1
		continue

elapsed = time.clock() - start
print "task done. time used: ", elapsed, "seconds"
## 释放内存
outputFile.close()
del resultAfterSort
exit()
