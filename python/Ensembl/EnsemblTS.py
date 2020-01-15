# coding=utf-8
# pzw
# 20200114

# database : ftp://ftp.ensembl.org/pub/grch37/current/fasta/homo_sapiens/cdna/Homo_sapiens.GRCh37.cdna.all.fa.gz

import sys
import gzip
import argparse

def ensemblTranscriptDB(db):
	grch37 = gzip.open(db, "rb")
	dbDict = {}
	for line in grch37.readlines():
		if line.startswith(">"):
			ts = line.split(">")[1].split(" ")[0]
			locations = line.split(" ")[2].split(":")
			location = locations[2] + ":" + locations[3] + "-" + locations[4]

			dbDict[ts] = [location, ""]

		else:
			dbDict[ts][1] += line.strip()
	grch37.close()
	return dbDict


def ensemblTranscriptSearch(transcriptID, dbDict):
	return dbDict[transcriptID]

def main(transcriptFile, outputFile, transcriptID, database):
	TsDB = ensemblTranscriptDB(database)
	
	if transcriptID:
		try:
			print ensemblTranscriptSearch(transcriptID, TsDB)
		except:
			print "转录本不存在"

	if transcriptFile:
		inputFile = open(transcriptFile, "r")
		output = open(outputFile, "w")
		for line in inputFile:
			ID = line.strip()
			try:
				outputResults = ensemblTranscriptSearch(ID, TsDB)
			except:
				print "未找到转录本： " + ID
				print "尝试寻找其他版本"
				for i in TsDB.keys():
					if i.split(".")[0] in ID:
						print "查找到相近版本： " + i + " 已写入，请手动调整"
						outputResults = ensemblTranscriptSearch(i, TsDB)

			output.write("\t".join(outputResults) + "\t" + ID + "\n")

		output.close()
		inputFile.close()

if __name__ == "__main__":
	parser = argparse.ArgumentParser(description="ensembl transcriptID search",
		prog="EnsemblTS.py",
		usage="python EnsemblTS.py -i <inputFile> -o <outputFile>")
	group = parser.add_mutually_exclusive_group()
	parser.add_argument("-v", "--version", action="version",
		version="Version 0.1 20200115")
	group.add_argument("-i", "--input", type=str,
		help="Input the file with one ensembl transcript ID on each line")
	parser.add_argument("-o", "--output", type=str,
		help="output file")
	group.add_argument("-I", "--ID", type=str,
		help="input a ensembl transcript ID. Example:ENST00000360863.6")
	parser.add_argument("-d", "--database", type=str,
		help="databases, default=hg19 cdna",
		default="Homo_sapiens.GRCh37.cdna.all.fa.gz")
	if len(sys.argv[1:]) == 0:
		parser.print_help()
		parser.exit()
	args = parser.parse_args()
	if args.input and args.ID:
		sys.exit("-i and -I must be only one!")
	main(transcriptFile=args.input, outputFile=args.output, transcriptID=args.ID, database=args.database)
