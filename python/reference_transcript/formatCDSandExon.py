# coding=utf-8
# pzw
# 20220623
# 规则化CDS及exon bed
# bed文件来源于UCSC

import sys

bedFile = sys.argv[1]
outputFile = sys.argv[2]

# databases
refTranscript = open("refTranscript.txt", "r")
refDict = {}
for r in refTranscript:
    if not r.startswith("#"):
        rs = r.split("\t")
        transcript = rs[1].split(".")[0]
        refDict[transcript] = rs[0]
refTranscript.close()

# bed file
bed = open(bedFile, "r")
output = open(outputFile, "w")
for line in bed:
    lines = line.split("\t")
    chrom = lines[0]
    start = lines[1]
    end = lines[2]
    name = lines[3]
    score = lines[4]
    strand = lines[5]
    if "_" in chrom:
        continue
    
    transcript = name.split(".")[0]
    if not transcript in refDict.keys():
        continue
    else:
        name = refDict[transcript]
    
    output.write("\t".join([chrom, start, end, name, transcript, strand]))

output.close()
bed.close()


