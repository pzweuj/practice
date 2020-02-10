# coding=utf-8
# pzw
# 20200210

import os
import sys
import argparse
import pandas as pd

# 多序列比对
def clustal(fasta, output, threads):
	clustalo = "/software/clustalo-1.2.4/clustalo"
	cmd = """
		{clustalo} -i {fasta} -o {output} --threads {threads}
	""".format(clustalo=clustalo, fasta=fasta, output=output, threads=threads)
	os.system(cmd)

# 结果整理为字典
def fastaDict(results):
	resultsFile = open(results, "r")
	resultsDict = {}
	for line in resultsFile:
		if line.startswith(">"):
			name = line.replace(">", "").split("\n")[0]
			resultsDict[name] = ""
		else:
			resultsDict[name] += line.split("\n")[0]
	resultsFile.close()

	return resultsDict

# 结果输出为文本
def rowText(results, outputFiles):
	output = open(outputFiles, "w")

	SeqDict = fastaDict(results)
	for i in SeqDict.keys():
		seq = i + "\t" + SeqDict[i] + "\n"
		output.write(seq)

	output.close()

# 结果输出到excel表格
def excelType(results, xlsxFile):

	xlsxDict = {}
	SeqDict = fastaDict(results)

	length = len(SeqDict.values()[0])
	xlsxDict["Index"] = []

	for i in range(length):
		indexNum = i + 1
		xlsxDict["Index"].append(indexNum)

	for j in SeqDict.keys():
		xlsxDict[j] = []
		for s in SeqDict[j]:
			xlsxDict[j].append(s)

	
	df = pd.DataFrame(xlsxDict)
	order = df.columns.values.tolist()
	order.remove("Index")
	order.insert(0, "Index")
	order.remove("MN908947")
	order.insert(1, "MN908947")
	df = df[order]
	df.to_excel(xlsxFile, index=False)
	print "done"

# 获得序列位置
def getIndex(results, primer):
	SeqDict = fastaDict(results)
	start = len(SeqDict["MN908947"].split(primer)[0]) + 1
	end = start + len(primer) - 1
	s = str(start) + "\t" + str(end)
	return s

# 大写
def seq_upper(seq):
	return seq.upper()

# 反向
def reverse(seq):
	return seq[::-1]

# 互补
def complement(seq):
	seqUpper = seq_upper(seq)
	seq_comp = seqUpper.replace("A", "t").replace("G", "c").replace("T", "a").replace("C", "g")
	seq_comp_degen = seq_comp.replace("R", "y").replace("Y", "r").replace("M", "k").replace("K", "m")
	seq_comp_degen2 = seq_comp_degen.replace("H", "d").replace("D", "h").replace("B", "v").replace("V", "b")
	return seq_upper(seq_comp_degen2)

# 反向互补
def complementReverse(seq):
	revSeq = reverse(seq)
	comp_Seq = complement(revSeq)
	return comp_Seq



def main(fasta, results, threads, output, primer, method):
	if method == "clustal":
		clustal(fasta, results, threads)
	elif method == "text":
		rowText(results, output)
	elif method == "xlsx":
		excelType(results, output)
	elif method == "primer":
		print getIndex(results, primer)
	elif method == "primer_com":
		print complementReverse(primer)
	else:
		print "please input a method: clustal/text/xlsx/primer/primer_com"

if __name__ == "__main__":
	parser = argparse.ArgumentParser(
		description="2019-nCoV analysis",
		prog="2019-nCoV.py",
		usage="""
			# 多序列比对
			python wuhan_2019-nCoV.py -f <all.fasta> -r <clustal.fasta> -t <threads> -m clustal
			# 文本格式结果输出
			python wuhan_2019-nCoV.py -r <clustal.fasta> -o <output.txt> -m text
			# excel格式结果输出
			python wuhan_2019-nCoV.py -r <clustal.fasta> -o <output.xlsx> -m xlsx
			# 引物位置检索
			python wuhan_2019-nCoV.py -r <clustal.fasta> -p <primer> -m primer
			# 引物反向互补
			python wuhan_2019-nCoV.py -p <primer> -m primer_com
		""",
		formatter_class=argparse.RawTextHelpFormatter
	)
	parser.add_argument("-v", "--version", action="version",
		version="Version 0.1 20200210")
	parser.add_argument("-f", "--fasta", type=str,
		help="在 https://bigd.big.ac.cn/ncov/genome/sequence/download/all 下载的all.fasta文件",
		default="")
	parser.add_argument("-r", "--results", type=str,
		help="多序列比对输出结果文件", default="")
	parser.add_argument("-t", "--threads", type=str,
		help="多序列比对使用线程数", default="8")
	parser.add_argument("-o", "--output", type=str,
		help="输出的文本结果或excel表格", default="")
	parser.add_argument("-p", "--primer", type=str,
		help="引物或探针序列", default="")
	parser.add_argument("-m", "--method", type=str,
		help="选择方法，可选clustal/text/xlsx/primer/primer_com", default="")
	if len(sys.argv[1:]) == 0:
		parser.print_help()
		parser.exit()
	args = parser.parse_args()
	main(fasta=args.fasta, results=args.results, threads=args.threads, output=args.output, primer=args.primer, method=args.method)