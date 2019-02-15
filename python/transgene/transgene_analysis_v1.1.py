# encoding:utf-8
# pzw
# 20190214
#####
# update log
# v0.2 更新了soft clip与map 的合并方式，根据先后顺序输出序列
# v0.3 更新阈值设置
# v0.4 加入正负链提示
# v1.0 读入文件为bam，自动读取质粒名称，自动判断是否比对到反义链
# v1.1 改成传参模式，统计总reads数的方法调整
#####
import re
import bamnostic
import time
import sys

#################################
bamFile = bamnostic.AlignmentFile(sys.argv[1], "rb")
results = open(sys.argv[2], "w")
basesCover = int(sys.argv[3])
showBases = int(sys.argv[4])			# 当设置为0时显示所有碱基；最终结果根据完全相同的值进行合并，设置值越大（0除外），结果越多
#################################


## 提示信息
print "bases Cover set: " + str(basesCover)
if showBases == 0:
	print "show all bases"
else:
	print "show " + str(showBases) + " bases"
print "-----------------------------------"
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
	readSplit = str(read).split("\t")
	hostChr = readSplit[2]
	readName = read.read_name
	cigar = read.cigarstring
	hostCigarSplit = re.findall(r"[0-9]+|[a-z]+", cigar)
	tagSA = str(read).split("SA:Z:")[1].split(",")
	insertCigarSplit = re.findall(r"[0-9]+|[a-z]+", tagSA[3])

	readCount.add(read.read_name)

	## 过滤以质粒为主要比对的reads
	if hostChr == plasmidName:
		print "remove read: " + readName + " of " + plasmidName
		continue

	## 过滤比对位置多于2的reads
	elif len(hostCigarSplit) > 2:
		print "remove read " + readName +": map to more than 2 locations"
		continue

	## 过滤覆盖区域未达设定碱基数的reads
	elif int(hostCigarSplit[0]) < basesCover or int(hostCigarSplit[1]) < basesCover:
		print "remove read " + readName + ": filter by basesCover"
		continue

	## 过滤输出bam文件时错误输出的reads
	elif tagSA[0] != plasmidName:
		print "remove read " + readName + ": wrong output"
		continue

	## 解析过滤后的reads
	else:
		print "analyse " + readName
		readName_count += 1
		seq = read.query_sequence
		hostSeq = read.query_alignment_sequence

		## 利用flag值判断正负链
		flags = bamnostic.utils.flag_decode(read.flag)
		if str(flags).__contains__("(16, 'read reverse strand')"):
			hostStrand = "-"
		else:
			hostStrand = "+"

		## 获得嵌合部分的位置与方向
		insertPos_tmp = tagSA[1]
		insertStrand = tagSA[2]
		
		## 判断宿主序列与质粒序列的合成方式
		if seq.split(hostSeq)[0] == "":
			insertSeq = seq.split(hostSeq)[1]
			interType = hostChr + "/" + plasmidName
			interStrand = hostStrand + "/" + insertStrand
			hostPos = int(readSplit[3]) + len(hostSeq) - 1
			insertPos = int(insertPos_tmp) + len(hostSeq) - int(insertCigarSplit[0])
			if showBases == 0:
				junctionSeq = hostSeq + "@" + insertSeq.lower()
			else:
				junctionSeq = hostSeq[-showBases:] + "@" + insertSeq[0:showBases].lower()
		else:
			insertSeq = seq.split(hostSeq)[0]
			interType = plasmidName + "/" + hostChr
			interStrand = insertStrand + "/" + hostStrand
			hostPos = int(readSplit[3])
			insertPos = int(insertPos_tmp) + len(insertSeq) + len(hostSeq) - int(insertCigarSplit[1]) - 1
			if showBases == 0:
				junctionSeq = insertSeq.lower() + "@" + hostSeq
			else:
				junctionSeq = insertSeq[-showBases:].lower() + "@" + hostSeq[0:showBases]

		## 输出
		outputList = [plasmidName, str(insertPos), hostChr, str(hostPos), interType, interStrand, junctionSeq]
		output = "\t".join(outputList)
		outputCount.append(output)

## 将相同的输出合并，统计总数
outputDict = {}
for output_tmp in outputCount:
	outputDict[output_tmp] = outputCount.count(output_tmp)

## 输出结果到文件
results.write("plasmidName\tplasmidPos\thostChr\thostPos\tintegrationType\tinterStrand\tfusionJunctionSeq\tcounts\n")
for output_final in outputDict.keys():
	result = output_final + "\t" + str(outputDict[output_final]) + "\n"
	results.write(result)

results.close()
bamFile.close()
elapsed = time.clock() - start

## 输出最终过滤结果
print "------------------------------------------"
print "input reads: " + str(len(readCount))
print "output reads: " + str(readName_count)
print str("%.2f" % (float(readName_count) / float(len(readCount)) * 100.00)) + "% reads pass"
print "task done. time used: ", elapsed, "seconds"
exit()