# coding=utf-8
# pzw
# 20200319

import os
import sys
import argparse

def GetVariantMatrix(AnnovarFilePath):
	file_list = os.listdir(AnnovarFilePath)
	changeList = []

	# 先遍历所有结果，获得所有变异位点
	for files in file_list:
		print(files)
		f = open(AnnovarFilePath + "/" + files, "r")
		for line in f:
			if line.startswith("Chr\tStart"):
				continue
			else:
				l = line.split("\t")
				chrom = l[0]
				start = l[1]
				end = l[2]
				ref = l[3]
				alt = l[4]
				exonicOrNot = l[5]
				gene = l[6]
				nonsynon = l[8]
				changeInfo = l[9]
				clinvar = l[48]

				# 筛选阈值
				bamDepth = l[74]
				if int(bamDepth) <= 500:
					continue

				# 甲状腺癌中，TERT为启动子区域突变，较为特殊
				if gene == "TERT":
					loc = start[4:]
					output = "\t".join([gene, chrom + ":" + start + "-" + end, "-", exonicOrNot, loc + ref + ">" + alt, "-", clinvar])
					changeList.append(output)

				# 其他突变只保留非同义突变
				else:
					if ("exonic" in exonicOrNot) or ("splicing" in exonicOrNot):
						if ("nonsynonymous" in nonsynon) or ("frameshift" in nonsynon):
							cc = changeInfo.split(":")
							transcript = cc[1]
							exon = cc[2]
							baseChange = cc[3]
							if len(cc) > 4:
								aachange = cc[4]
							else:
								aachange = "-"

							if "," in aachange:
								aachange = aachange.split(",")[0]

							output = "\t".join([gene, chrom + ":" + start + "-" + end, transcript, exon, baseChange, aachange, clinvar])
							changeList.append(output)
					else:
						continue
		f.close()
	# 获得的变异位点去重并排序
	changeList = sorted(list(set(changeList)))
	return changeList

def MakeOutputMatrix(AnnovarFilePath, outputFile):
	outputFile = open(outputFile, "w")
	# 输出的头信息
	line0 = "changeID"
	line1 = "Gene"
	line2 = "Location"
	line3 = "transcript"
	line4 = "exon"
	line5 = "baseChange"
	line6 = "aaChange"
	line7 = "clinvar"

	# 统计获得变异位点数并生成唯一编号
	changeList = GetVariantMatrix(AnnovarFilePath)
	n = 0
	for i in changeList:
		n += 1
		line0 = line0 + "\t" + str(n)
		ii = i.split("\t")
		line1 = line1 + "\t" + ii[0]
		line2 = line2 + "\t" + ii[1]
		line3 = line3 + "\t" + ii[2]
		line4 = line4 + "\t" + ii[3]
		line5 = line5 + "\t" + ii[4]
		line6 = line6 + "\t" + ii[5]
		line7 = line7 + "\t" + ii[6]


	outputFile.write(line0 + "\n")
	outputFile.write(line1 + "\n")
	outputFile.write(line2 + "\n")
	outputFile.write(line3 + "\n")
	outputFile.write(line4 + "\n")
	outputFile.write(line5 + "\n")
	outputFile.write(line6 + "\n")
	outputFile.write(line7 + "\n")

	# 第二次遍历文件
	file_list = os.listdir(AnnovarFilePath)
	amount = len(file_list)
	file_n = 0
	for files in file_list:
		f = open(AnnovarFilePath + "/" + files, "r")
		file_n += 1
		stringO = "正在分析第" + str(file_n) + "个样本(总计" + str(amount) + "个): " + files
		print(stringO)

		Results_f = []
		# 遍历变异位点，逐一填写结果
		for insertInfo in changeList:
			insertX = ""
			fx = open(AnnovarFilePath + "/" + files, "r")
			for line in fx:
				if line.startswith("Chr\tStart"):
					continue
				else:
					l = line.split("\t")
					chrom = l[0]
					start = l[1]
					end = l[2]
					nonsynon = l[8]
					changeInfo = l[9]
					bamDepth = l[74]

					if int(bamDepth) <= 500:
						continue

					otherInfo = l[84]
					AF = otherInfo.split(":")[2]
					if "," in AF:
						AF = AF.split(",")[1]
					DP = bamDepth
					MAF = "%.4f" % (float(AF) / int(DP))

					checkPoint = chrom + ":" + start + "-" + end

					# 检查是否存在
					if checkPoint == insertInfo.split("\t")[1]:

						if len(changeInfo.split(":")) > 4:
							aachange = changeInfo.split(":")[4]
						else:
							aachange = "-"

						if "," in aachange:
							aachange = aachange.split(",")[0]

						if aachange != insertInfo.split("\t")[5]:
							continue

						# 检查是否满足要求，AF>10，DP>=5000，MAF>0.2%
						if int(AF) >= 5:
							if float(MAF) >= 0.001:
								insertX = MAF + "(" + str(AF) + "/" + DP + ")"
			fx.close()
			Results_f.append(insertX)
		outputFile.write(files.split(".")[0] + "\t" + "\t".join(Results_f) + "\n")

		f.close()
	outputFile.close()
	print("Task Done")

def main(AnnovarFilePath, outputFile):
	GetVariantMatrix(AnnovarFilePath)
	MakeOutputMatrix(AnnovarFilePath, outputFile)

if __name__ == "__main__":
	parser = argparse.ArgumentParser(description="Thyroid nodule variants matrix",
		prog="thyroid_matrix.py",
		usage="python3 thyroid_matrix.py -i <annovar reults path> -o <matrix file>")
	group = parser.add_mutually_exclusive_group()
	parser.add_argument("-v", "--version", action="version",
		version="Version 0.1 20200330")
	parser.add_argument("-i", "--input", type=str,
		help="Input the directory of the multianno file")
	parser.add_argument("-o", "--output", type=str,
		help="output the matrix file")
	if len(sys.argv[1:]) == 0:
		parser.print_help()
		parser.exit()
	args = parser.parse_args()
	main(AnnovarFilePath=args.input, outputFile=args.output)