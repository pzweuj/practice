# encoding:utf-8
# pzw
# 20190225
#####
# update log
# v2.1 方法调整
# v2.0 方法调整
# v1.6 更新统计方法，提升性能
# v1.5 使用内置函数
# v1.4 增加连接位点之间有两边均比对不上的碱基的过滤阈值，人为规定为10，目前未作为可选参数
# v1.3 修正过滤逻辑
# v1.2 部分连接位点之间有两边均比对不上的碱基，修复这部分的内容，调整结果呈现方式
# v1.1 改成传参模式，统计总reads数的方法调整
# v1.0 读入文件为bam，自动读取质粒名称，自动判断是否比对到反义链
# v0.4 加入正负链提示
# v0.3 更新阈值设置
# v0.2 更新了soft clip与map 的合并方式，根据先后顺序输出序列
#####
import re
import bamnostic
import time
import sys
from collections import Counter

#################################
bamFile = bamnostic.AlignmentFile(sys.argv[1], "rb")
results = open(sys.argv[2], "w")
basesCover = int(sys.argv[3])
showBases = int(sys.argv[4])		# 当设置为0时显示所有碱基；最终结果根据完全相同的值进行合并，设置值越大（0除外），结果越多
maxGap = 10
#################################

##########  test  ###############
# bamFile = bamnostic.AlignmentFile("H04730D.CED1484.final.bam", "rb")
# results = open("results.txt", "w")
# basesCover = 46
# showBases = 10
# maxGap = 10
#################################

## 提示信息
print("bases Cover set: " + str(basesCover))
if showBases == 0:
	print("show all bases")
else:
	print("show " + str(showBases) + " bases")
print("-----------------------------------")
start = time.clock()

## 读取头信息
head = bamFile.header
plasmidName = head["SQ"][-1]["SN"]

outputCount = []

## reads统计
readCount = set()
readName_count = 0

for read in bamFile:
	## 预处理
	readName = read.read_name
	mateName = read.next_reference_name
	cigar = read.cigarstring
	insertCigarNum = re.findall(r"[0-9]+|[a-z]+", cigar)
	readCount.add(read.read_name)
	if str(read).__contains__("SA:Z"):
		chimera = "positive"
	else:
		chimera = "negative"

	## 过滤比对到多个位置的reads
	if len(insertCigarNum) > 2:
		print "remove read " + readName +": map to more than 2 locations"
		continue

	## 分析
	if mateName == plasmidName:
		if chimera == "positive":
			tagSASplit = read.get_tag("SA").split(",")
			hostChr = tagSASplit[0]
			hostPos_tmp = tagSASplit[1]
			hostCigar = tagSASplit[3]
			hostCigarNum = re.findall(r"[0-9]+|[a-z]+", hostCigar)

			## 过滤嵌合部分比对到质粒的reads
			if hostChr == plasmidName:
				print("remove read " + readName +": both map to plasmid")
				continue
		else:
			continue	
	else:
		continue

	## 过滤比对到多个位置的reads
	if len(hostCigarNum) > 2:
		print("remove read " + readName +": map to more than 2 locations")
		continue

	## 宿主序列的合成方式
	hostCigarDirect = hostCigar.split("M")
	if hostCigarDirect[1] == "":
		hostMapAt = "right"
		hostSoft = hostCigar.split("S")[0]
		hostMap = hostCigar.split("S")[1].split("M")[0]
	else:
		hostMapAt = "left"
		hostSoft = hostCigarDirect[1].split("S")[0]
		hostMap = hostCigarDirect[0]

	## 质粒序列的合成方式
	insertCigarDirect = cigar.split("M")
	if insertCigarDirect[1] == "":
		insertMapAt = "right"
		insertSoft = insertCigarDirect[0].split("S")[0]
		insertMap = insertCigarDirect[0].split("S")[1]
	else:
		insertMapAt = "left"
		insertSoft = insertCigarDirect[1].split("S")[0]
		insertMap = insertCigarDirect[0]

	## 过滤宿主序列与质粒序列属于包含关系的reads
	if hostMapAt == insertMapAt:
		print("remove read " + readName + ": filter by inclusion")
		continue

	## 过滤覆盖区域未达设定碱基数的reads
	if int(hostMap) < 20 or int(insertMap) < basesCover:
		print("remove read " + readName + ": filter by basesCover")
		continue

	## 过滤掉中间间隔太大的reads
	if (int(insertMap) + maxGap) < int(hostSoft):
		print("remove read " + readName + ": gap too large")
		continue

	## 解析过滤后的reads
	else:
		print("analyse " + readName)
		readName_count += 1
		seq = read.query_sequence
		insertSeq = read.query_alignment_sequence

		## 利用flag值判断正负链
		# if read.is_reverse:
		# 	insertStrand = "-"
		# else:
		# 	insertStrand = "+"

		if insertMapAt == "left":
			hostSeq = seq[-int(hostMap):]
			interType = plasmidName + "/" + hostChr
			insertPos = read.query_alignment_start + read.query_alignment_length
			hostPos = int(hostPos_tmp)
			if showBases == 0:
				junctionSeq = insertSeq.lower() + "@" + hostSeq
			else:
				junctionSeq = insertSeq[-showBases:].lower() + "@" + hostSeq[0:showBases]
		else:
			hostSeq = seq[0:int(hostMap)]
			interType = hostChr + "/" + plasmidName
			insertPos = read.query_alignment_start + 1
			hostPos = int(hostPos_tmp) + int(hostMap) - 1
			if showBases == 0:
				junctionSeq = hostSeq + "@" + insertSeq.lower()
			else:
				junctionSeq = hostSeq[-showBases:] + "@" + insertSeq[0:showBases].lower()

		## 输出
		outputList = [
			plasmidName,
			str(insertPos),
			hostChr,
			str(hostPos),
			interType,
			junctionSeq
		]
		output = "\t".join(outputList)
		outputCount.append(output)

# 将相同的输出合并，统计总数
outputDict = Counter(outputCount)
del outputCount

## 输出结果到文件
results.write("plasmidName\tplasmidPos\thostChr\thostPos\tintegrationType\tfusionJunctionSeq\tcounts\n")
for output_final in outputDict.keys():
	result = output_final + "\t" + str(outputDict[output_final]) + "\n"
	results.write(result)

results.close()
bamFile.close()
elapsed = time.clock() - start

## 输出最终过滤结果
print("------------------------------------------")
print("input reads: " + str(len(readCount)))
print("output reads: " + str(readName_count))
print(str("%.2f" % (float(readName_count) / float(len(readCount)) * 100.00)) + "% reads pass")
print("task done. time used: ", elapsed, "seconds")
exit()
