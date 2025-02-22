# encoding:utf-8
# pzw
# 20190214
#####
# update log
# v0.2 更新了soft clip与map 的合并方式，根据先后顺序输出序列
# v0.3 更新阈值设置
# v0.4 加入正负链提示
#####
import re

## 输入sam文件
samFile = open("transgene.final.R.sam", "r")
## 质粒名称
plasmidName = "pGRG36"
## 输出结果
results = open("results.R.txt", "w")
## 设定输出序列的长度，若无输出，需要降低设置值
seqLength = 20
## 阈值
basesCover = 20
## 方向
strand = "-"


outputCount = []
for line in samFile:
	if line.startswith("@"):
		continue
	else:
		lineAfterSplit = line.split("\t")
		splitReadsLength = re.findall(r"[0-9]+|[a-z]+", lineAfterSplit[5])

		seq = lineAfterSplit[9]
		hostStrand = strand

		if len(splitReadsLength) == 2:
			readName = lineAfterSplit[0]
			hostName = lineAfterSplit[2]
			hostStart = lineAfterSplit[3]
			
			## 避免不同的结果中SA tag的位置不同
			for tag in lineAfterSplit:
				if tag.__contains__("SA:Z"):
					insertName = tag.split(":")[2].split(",")[0]
					insertStart = tag.split(":")[2].split(",")[1]
					insertStrand = tag.split(":")[2].split(",")[2]

					mapLengthSplit = lineAfterSplit[5].split("M")
					if mapLengthSplit[1] == "":
						mapLength = mapLengthSplit[0].split("S")[1]
						softClipLength = mapLengthSplit[0].split("S")[0]
						mapSeq = seq[int(softClipLength):]
						softSeq = seq[0:int(softClipLength)]
						interType = insertName + "/" + hostName
						interStrand = insertStrand + "/" + hostStrand
						if hostName == plasmidName:
							fusionJunctionSeq = softSeq[0-int(seqLength):] + "@" + mapSeq[0:seqLength].lower()
						else:
							fusionJunctionSeq = softSeq[0-int(seqLength):].lower() + "@" + mapSeq[0:seqLength]
					else:
						mapLength = mapLengthSplit[0]
						softClipLength = mapLengthSplit[1].split("S")[0]
						mapSeq = seq[0:int(mapLength)]
						softSeq = seq[int(mapLength):]
						interType = hostName + "/" + insertName
						interStrand = hostStrand + "/" + insertStrand
						if hostName == plasmidName:
							fusionJunctionSeq = mapSeq[0-int(seqLength):].lower() + "@" + softSeq[0:seqLength]
						else:
							fusionJunctionSeq = mapSeq[0-int(seqLength):] + "@" + softSeq[0:seqLength].lower()

					hostEnd = int(hostStart) + int(mapLength) - 1

					## 避免输出染色体自身的junction reads
					if hostName == insertName:
						continue

					## 避免多比对情况出现
					elif hostName != plasmidName and insertName != plasmidName:
						continue

					## 使用阈值去除未达标的reads
					elif mapLength < basesCover or softClipLength < basesCover:
						continue
					else:
						insertEnd = int(softClipLength) + int(insertStart) - 1

						if hostName == plasmidName:
							plasmidPos = hostEnd
							hostChr = insertName
							hostPos = insertStart
						else:
							plasmidPos = insertStart
							hostChr = hostName
							hostPos = hostEnd

						outputList = [plasmidName, str(plasmidPos), hostChr, str(hostPos), interType, interStrand, fusionJunctionSeq]
						output = "\t".join(outputList)

						outputCount.append(output)

outputDict = {}
for i in outputCount:
	outputDict[i] = outputCount.count(i)

results.write("plasmidName\tplasmidPos\thostChr\thostPos\tintegrationType\tinterStrand\tfusionJunctionSeq\tcounts\n")
for j in outputDict.keys():
	result = j + "\t" + str(outputDict[j]) + "\n"
	results.write(result)

results.close()
samFile.close()