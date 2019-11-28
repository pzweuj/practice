# coding=utf-8
# pzw
# 20190627
# v0.2

import sys
import argparse


# 获得参考序列
def getRefSeq(fasta):
	ref = open(fasta, "r")
	refSeq = ref.readlines()[1]
	ref.close()
	return refSeq

# 导入vcf，获得突变序列
def BaseChange(vcfFile):
	vcf = open(vcfFile, "r")
	baseChange = {}
	for line in vcf:
		if line.startswith("#"):
			continue
		else:
			lineAfterSplit = line.split("\t")
			pos = lineAfterSplit[1]
			ref = lineAfterSplit[3]
			alt = lineAfterSplit[4]
			DP = lineAfterSplit[7].split("DP=")[1].split(";")[0]
			baseChange[pos] = [ref, alt, DP]

	vcf.close()
	return baseChange

# 翻译
def translate(seq):
	table = { 
		"ATA":"I", "ATC":"I", "ATT":"I", "ATG":"M",
		"ACA":"T", "ACC":"T", "ACG":"T", "ACT":"T",
		"AAC":"N", "AAT":"N", "AAA":"K", "AAG":"K",
		"AGC":"S", "AGT":"S", "AGA":"R", "AGG":"R",
		"CTA":"L", "CTC":"L", "CTG":"L", "CTT":"L",
		"CCA":"P", "CCC":"P", "CCG":"P", "CCT":"P",
		"CAC":"H", "CAT":"H", "CAA":"Q", "CAG":"Q",
		"CGA":"R", "CGC":"R", "CGG":"R", "CGT":"R",
		"GTA":"V", "GTC":"V", "GTG":"V", "GTT":"V",
		"GCA":"A", "GCC":"A", "GCG":"A", "GCT":"A",
		"GAC":"D", "GAT":"D", "GAA":"E", "GAG":"E",
		"GGA":"G", "GGC":"G", "GGG":"G", "GGT":"G",
		"TCA":"S", "TCC":"S", "TCG":"S", "TCT":"S",
		"TTC":"F", "TTT":"F", "TTA":"L", "TTG":"L",
		"TAC":"Y", "TAT":"Y", "TAA":"_", "TAG":"_",
		"TGC":"C", "TGT":"C", "TGA":"_", "TGG":"W",
	}
	protein = ""
	if len(seq) % 3 == 0:
		for i in range(0, len(seq), 3):
			codon = seq[i: i + 3]
			protein += table[codon]
	else:
		print ("Warning, sequence can not be translated completely")
		for i in range(0, len(seq), 3):
			if i + 3 < len(seq):
				codon = seq[i: i + 3]
				protein += table[codon]
	return protein

# 获得变异后氨基酸序列
def ChangeAnimo(vcfFile, reference, depthCutOff, skipIns):

	baseSeq = []
	for i in getRefSeq(reference):
		baseSeq.append(i)

	baseChange = BaseChange(vcfFile)

	for p in baseChange.keys():
		if int(baseChange[p][1]) >= depthCutOff:
			if skipIns:
				if len(baseChange[p][1]) > 1:
					continue
				elif len(baseChange[p][0]) > 1:
					continue
				else:
					baseSeq[int(p)-1] = baseChange[p][1]
			else:
				baseSeq[int(p)-1] = baseChange[p][1]

	changeSeq = "".join(baseSeq)
	changeAnimo = translate(changeSeq)

	return changeSeq, changeAnimo

# 获得变异氨基酸位点
def changeLocation(vcfFile, fasta, depthCutOff, skipIns):

	seq = translate(getRefSeq(fasta))
	changeSeq = ChangeAnimo(vcfFile, fasta, depthCutOff, skipIns)[1]

	changeLocationDict = {}
	for i in range(len(seq)):
		if seq[i] == changeSeq[i]:
			continue
		else:
			changeLocationDict[str(i+1)] = [seq[i], changeSeq[i]]
	return changeLocationDict

# 获得表格结果和药物结果
def drugLocation(vcfFile, fasta, depthCutOff, drugDB, skipIns):
	drug = open(drugDB, "r")
	changLoc = changeLocation(vcfFile, fasta, depthCutOff, skipIns)

	drugList = []
	outputList = []

	for line in drug:
		if line.startswith("#"):
			continue
		else:
			lines = line.split("\t")
			loca = lines[0]
			ref = lines[1]
			# alt = lines[2]
			drugs = lines[3].split("\n")[0]

			if loca in changLoc.keys():
				outputList.append("rt" + changLoc[loca][0] + loca + changLoc[loca][1])
				drugList.append(drugs)
			else:
				outputList.append("未突变")
	drug.close()

	drugList = list(set(drugList))

	# 恩替卡韦检验，当存在rtM204V以及rtL180M且有至少一个其他位点变异时耐药
	if "恩替卡韦" in drugList:
		if "rtM204V" not in outputList or "rtL180M" not in outputList:
			drugList.remove("恩替卡韦")
		else:
			drugETV = open(drugDB, "r")
			n = 0
			for lineETV in drugETV:
				if "恩替卡韦" in lineETV:
					if lineETV.split("\t")[0] in changLoc.keys():
						n += 1
			if n <= 2:
				drugList.remove("恩替卡韦")

	# 替比夫定检验，当rtM204V单独存在，不耐药
	if "替比夫定" in drugList:
		if "rtM204V" in outputList and "rtL180M" not in outputList:
			drugList.remove("替比夫定")

	return drugList, outputList
	

def main(inputVcf, reference, drugDB, depth, method, skipIns):
	drugFind = drugLocation(inputVcf, reference, depth, drugDB, skipIns)[0]
	changeList = drugLocation(inputVcf, reference, depth, drugDB, skipIns)[1]
	changeSeq = ChangeAnimo(inputVcf, reference, depth, skipIns)[0]
	
	if method == "drug":
		if len(drugFind) > 0:
			print "，".join(list(drugFind))
		else:
			print "未发现耐药突变"

	elif method == "list":
		for i in changeList:
			print i

	elif method == "seq":
		# S gene sequence
		print changeSeq[25:706]

	else:
		print "wrong method, please use 'drug' or 'list'."

if __name__ == "__main__":
	parser = argparse.ArgumentParser(description="HBV Drug @PZW",
		prog="HBV_drug.py",
		usage="python HBV_drug.py -i <inputVcf> -c <cutoff> -r <reference> -d <drug database> -m <method>")
	parser.add_argument("-v", "--version", action="version",
		version="Version 0.3 20190813")
	parser.add_argument("-i", "--input", type=str,
		help="Input the bcftools result file")
	parser.add_argument("-c", "--cutoff", type=int,
		help="min depth")
	parser.add_argument("-r", "--reference", type=str,
		help="reference file")
	parser.add_argument("-d", "--drug", type=str,
		help="drug database")
	parser.add_argument("-m", "--method", type=str,
		help="list: output change list, drug: output drug, seq: output RT region", default="drug")
	parser.add_argument("-s", "--skip", type=bool,
		help="skip insertion", default=False)
	if len(sys.argv[1:]) == 0:
		parser.print_help()
		parser.exit()
	args = parser.parse_args()
	main(inputVcf=args.input, reference=args.reference, drugDB=args.drug, depth=args.cutoff, method=args.method, skipIns=args.skip)
