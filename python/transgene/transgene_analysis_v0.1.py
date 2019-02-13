# encoding:utf-8
# pzw
# 20190213
import re

samFile = open("transgene.final.sam", "r")
plasmidName = "pGRG36"
results = open("results.txt", "w")

## 设定输出序列的长度，若无输出，需要降低设置值
seqLength = 20

outputCount = []
for line in samFile:
	if line.startswith("@"):
		continue
	else:
		lineAfterSplit = line.split("\t")
		splitReadsLength = re.findall(r"[0-9]+|[a-z]+", lineAfterSplit[5])
		if len(splitReadsLength) == 2:
			readName = lineAfterSplit[0]
			hostName = lineAfterSplit[2]
			hostStart = lineAfterSplit[3]
			hostEnd = int(splitReadsLength[0]) + int(hostStart) - 1
			
			## 避免不同的结果中SA tag的位置不同
			for tag in lineAfterSplit:
				if tag.__contains__("SA:Z"):
					insertName = tag.split(":")[2].split(",")[0]
					insertStart = tag.split(":")[2].split(",")[1]

					## 避免输出染色体自身的junction reads
					if hostName == insertName:
						continue

					## 避免多比对情况出现
					elif hostName != plasmidName and insertName != plasmidName:
						continue
					else:
						insertEnd = int(splitReadsLength[1]) + int(insertStart) - 1

						seq = lineAfterSplit[9]
						hostSeq = seq[int(splitReadsLength[0])-seqLength:int(splitReadsLength[0])]
						insertSeq = seq[int(splitReadsLength[0])+1:int(splitReadsLength[0])+seqLength]

						interType = hostName + "/" + insertName

						if hostName == plasmidName:
							plasmidPos = hostEnd
							hostChr = insertName
							hostPos = insertStart
							hostSeq = hostSeq.lower()
						else:
							plasmidPos = insertStart
							hostChr = hostName
							hostPos = hostEnd
							insertSeq = insertSeq.lower()

						fusionJunctionSeq = hostSeq + "@" + insertSeq

						outputList = [plasmidName, str(plasmidPos), hostChr, str(hostPos), interType, fusionJunctionSeq]
						output = "\t".join(outputList)

						outputCount.append(output)

outputDict = {}
for i in outputCount:
	outputDict[i] = outputCount.count(i)

results.write("plasmidName\tplasmidPos\thostChr\thostPos\tintegrationType\tfusionJunctionSeq\tcounts\n")
for j in outputDict.keys():
	result = j + "\t" + str(outputDict[j]) + "\n"
	results.write(result)

results.close()
samFile.close()