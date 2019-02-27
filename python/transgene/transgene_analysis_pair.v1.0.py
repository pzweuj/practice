# encoding:utf-8
# pzw
# 20190227
#####
# update log
# v1.0 传参模式
# v0.3 运行log，统计，补全注释
# v0.2 校正坐标
# v0.1
#####
import bamnostic
import time
import sys
import pandas as pd

#################################
bamFile = bamnostic.AlignmentFile(sys.argv[1], "rb")
results = sys.argv[2]
basesCover = int(sys.argv[3])
showBases = int(sys.argv[4])	# 当设置为0时显示所有碱基；最终结果根据完全相同的值进行合并，设置值越大（0除外），结果越多
hostBasesCover = 20

##########  test  ###############
# bamFile = bamnostic.AlignmentFile("H04730D.pair.final.bam", "rb")
# results = "results.paired.txt"
# basesCover = 40
# showBases = 10	# 当设置为0时显示所有碱基；最终结果根据完全相同的值进行合并，设置值越大（0除外），结果越多
# hostBasesCover = 20
#################################

### 反向互补 ###
# def DNA_complement(sequence):
# 	sequence = sequence.upper()
# 	sequence = sequence.replace("A", "t")
# 	sequence = sequence.replace("T", "a")
# 	sequence = sequence.replace("C", "g")
# 	sequence = sequence.replace("G", "c")
# 	return sequence.upper()
 
# def DNA_reverse(sequence):
# 	sequence = sequence.upper()
# 	return sequence[::-1]

## 提示信息
print "bases Cover set: " + str(basesCover)
if showBases == 0:
	print "show all bases"
else:
	print "show " + str(showBases) + " bases"
print "-----------------------------------"
start = time.clock()

read_count = 0

plasmidName = bamFile.header["SQ"][-1]["SN"]
df_host = pd.DataFrame(columns=["read_ID", "hostChr", "hostPos", "integrationType", "hostSeq"])
df_plasmid = pd.DataFrame(columns=["read_ID", "plasmidName", "insertPos", "insertSeq"])
for read in bamFile:

	## 读入基本信息
	read_ID = read.read_name
	chrName = str(read).split("\t")[2]
	nextName = read.next_reference_name
	read_count += 1

	## 当read比对到质粒
	if chrName == plasmidName:
		insertName = chrName
		hostName = nextName
		insertSeq = read.query_alignment_sequence

		## 过滤cover值
		if len(insertSeq) < basesCover:
			print read_ID, "filter by basesCover"
			continue

		## 判断是否是read1
		if read.is_read1:
			insertPos = read.query_alignment_start + read.query_alignment_length
			interType = insertName + "/" + hostName

			## 输出碱基数
			if showBases != 0:
				insertSeq = insertSeq[-showBases:]

		## read2
		else:
			insertPos = read.query_alignment_start + 1
			interType = hostName + "/" + insertName

			## 输出碱基数
			if showBases != 0:
				insertSeq = insertSeq[0:showBases]

		## 输出以质粒为主要比对的read到表格
		print "writing", read_ID
		df_plasmid.loc[len(df_plasmid.index)] = [read_ID, insertName, insertPos, insertSeq]
		plasmidReads_count = len(df_plasmid.index)

	## 当配对的read比对到质粒
	# 即该read比对到宿主，由于在生成bam文件时，筛选的是read1/2比对到不同染色体的read pair，因此不会存在read1/2同时比对到相同染色体的情况
	elif nextName == plasmidName:
		insertName = nextName
		hostName = chrName
		hostSeq = read.query_alignment_sequence

		## 宿主序列的阈值过滤
		if len(hostSeq) < hostBasesCover:
			print read_ID, "filter by basesCover"
			continue
		if read.is_read1:
			hostPos = read.query_alignment_start + read.query_alignment_length
			interType = hostName + "/" + insertName

			## 输出碱基数
			if showBases != 0:
				hostSeq = hostSeq[-showBases:]
		else:
			hostPos = read.query_alignment_start + 1
			interType = insertName + "/" + hostName

			## 输出碱基数
			if showBases != 0:
				hostSeq = hostSeq[0:showBases]

		## 输出以宿主为主要比对的read到表格
		print "writing", read_ID
		df_host.loc[len(df_host.index)] = [read_ID, hostName, hostPos, interType, hostSeq]
		hostReads_count = len(df_host.index)

	## 过滤read1/2均没有比对到质粒的
	else:
		print read_ID, "not interest"
		continue

bamFile.close()

## 合并pairend 数据
print "merging DataFrame"
df_final = pd.merge(df_host, df_plasmid, on="read_ID", how="inner")
df_final["fusionJunctionSeq"] = ""
del df_host
del df_plasmid
pairReads_count = len(df_final.index)

## 处理fusion Junction
print "merging sequence"
for i in df_final.index:
	if df_final.loc[i, "integrationType"].split("/")[0] == plasmidName:
		df_final.ix[i, "fusionJunctionSeq"] = df_final.loc[i, "insertSeq"].lower() + "@" + df_final.loc[i, "hostSeq"]
	else:
		df_final.ix[i, "fusionJunctionSeq"] = df_final.loc[i, "hostSeq"] + "@" + df_final.loc[i, "insertSeq"].lower()


df_final["counts"] = 1
## 下面这一行可完整输出所有的reads，需要时调整注释
# df_final = df_final[["read_ID", "plasmidName", "insertPos", "hostName", "hostPos", "integrationType", "fusionJunctionSeq"]]
df_final = df_final[["plasmidName", "insertPos", "hostChr", "hostPos", "integrationType", "fusionJunctionSeq", "counts"]]

## 去重统计
df_final = df_final.groupby(["plasmidName", "insertPos", "hostChr", "hostPos", "integrationType", "fusionJunctionSeq"]).sum().reset_index()
df_final.to_csv(results, sep="\t", index=False)
print "output done"
del df_final


## 输出统计值
print "read count: ", read_count
print "plasmid read count: ", plasmidReads_count
print "host read count: ", hostReads_count
print "read pairs count: ", pairReads_count

elapsed = time.clock() - start
print "task done. time used: ", elapsed, "seconds"
exit()