# coding=utf-8
# pzw
# 20190627
# v0.2

import sys
import argparse

def getRefSeq(fasta):
	ref = open(fasta, "r")
	refSeq = ref.readlines()[1]
	ref.close()
	return refSeq

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
			baseChange[pos] = [alt, DP]
	vcf.close()
	return baseChange


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

def ChangeAnimo(vcfFile, reference, depthCutOff):

	baseSeq = []
	for i in getRefSeq(reference):
		baseSeq.append(i)

	baseChange = BaseChange(vcfFile)

	for p in baseChange.keys():
		if int(baseChange[p][1]) >= depthCutOff:
			baseSeq[int(p)-1] = baseChange[p][0]

	changeSeq = "".join(baseSeq)
	changeAnimo = translate(changeSeq)

	return changeAnimo


def drugLocation(seq, drugDB):
	drug = open(drugDB, "r")

	res_drug = set()
	for line in drug:
		if line.startswith("#"):
			continue
		else:
			l = line.split("\t")
			loca = int(l[0])
			ref = l[1]
			alt = l[2]
			drugs = l[3].split("\n")[0]

			if seq[loca-1] == alt:
				res_drug.add(drugs)
				print (ref + str(loca) + alt), drugs

	drug.close()
	return res_drug

def main(inputVcf, reference, drugDB, depth):
	changeSeq = ChangeAnimo(inputVcf, reference, depth)
	drugFind = drugLocation(changeSeq, drugDB)
	
	if len(drugFind) > 0:
		print "耐药： ", ",".join(list(drugFind))
	else:
		print "未发现耐药突变"

if __name__ == "__main__":
	parser = argparse.ArgumentParser(description="HBV Drug @PZW",
		prog="HBV_drug.py",
		usage="python HBV_drug.py -i <inputVcf> -c <cutoff> -r <reference> -d <drug database>")
	parser.add_argument("-v", "--version", action="version",
		version="Version 0.2 20190627")
	parser.add_argument("-i", "--input", type=str,
		help="Input the bcftools result file")
	parser.add_argument("-c", "--cutoff", type=int,
		help="min depth")
	parser.add_argument("-r", "--reference", type=str,
		help="reference file")
	parser.add_argument("-d", "--drug", type=str,
		help="drug database")
	if len(sys.argv[1:]) == 0:
		parser.print_help()
		parser.exit()
	args = parser.parse_args()
	main(inputVcf=args.input, reference=args.reference, drugDB=args.drug, depth=args.cutoff)
