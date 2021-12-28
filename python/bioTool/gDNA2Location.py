# coding=utf-8
# pzw
# 20211228
# trans chrN:g.1000000A>T format to chrN    10000000    10000000    A   T

import os
import re

def gDNA2Loca(query):
    if ":" in query:
        querys = query.split(":")
        chrom = querys[0]
        locas = querys[1].replace("g.", "").replace("(", "").replace(")", "")

        locaSplit = locas.split("_")
        # 7578232A>T
        if len(locaSplit) == 1:
            if "del" in locas:
                changes = locaSplit[0].split("del")
                alt = "-"
                ref = changes[1]
                start = pos = changes[0]
                stop = str(int(pos) + len(alt) - 1)
            elif "dup" in locas:
                changes = locaSplit[0].split("dup")
                pos = changes[0]
                start = pos
                stop = "-"
                ref = "-"
                alt = changes[1]
            else:
                changes = locaSplit[0].split(">")
                alt = changes[1]
                pos_ref = re.findall(r'[A-Za-z]+|\d+', changes[0])
                pos = pos_ref[0]
                ref = pos_ref[1]
                start = pos
                stop = str(int(pos) + len(alt) - 1)
        elif len(locaSplit) == 2:
            if "ins" in locaSplit[1] and "del" in locaSplit[1]:
                pos = start = locaSplit[0]
                stop_change = locaSplit[1].split("del")
                change = stop_change[1]
                stop = stop_change[0]
                del_ins = change.split("ins")
                ref = del_ins[0]
                alt = del_ins[1]
            elif "ins" in locaSplit[1]:
                pos = start = locaSplit[0]
                stop_change = locaSplit[1].split("ins")
                stop = "-"
                alt = change = stop_change[1]
                ref = "-"
            elif "del" in locaSplit[1]:
                pos = start = locaSplit[0]
                stop_change = locaSplit[1].split("del")
                ref = stop_change[1]
                stop = stop_change[0]
                alt = "-"
            elif "dup" in locaSplit[1]:
                pos = locaSplit[0]
                stop_change = locaSplit[1].split("dup")
                ref = "-"
                alt = stop_change[1]
                start = stop_change[0]
                stop = "-"
            else:
                pos = start = locaSplit[0]
                stop = locaSplit[1]
                ref = ""
                alt = "-"

        else:
            chrom = start = stop = ref = alt = ""
    else:
        chrom = start = stop = ref = alt = ""

    if ref.isdigit():
        ref = ""

    return [chrom, start, stop, ref, alt]

def fixReference(chrom, start, stop, reference):
    cmd = """
        samtools faidx {reference} {chrom}:{start}-{stop}
    """.format(reference=reference, chrom=chrom, start=start, stop=stop)
    
    r = os.popen(cmd)
    results = r.readlines()
    seq = ""
    for line in results:
        if not line.startswith(">"):
            seq = seq + line.replace("\n", "")
    return seq

def fixFormat(query):
    reference = "/slurm/databases/hg19/ucsc.hg19.fasta"
    locaList = gDNA2Loca(query)
    fixChrom = locaList[0]
    print(locaList)
    if locaList[2] == "-":
        fixStart = locaList[1]
        fixStop = fixStart
        fixRef = fixReference(fixChrom, fixStart, fixStop, reference)
        fixAlt = fixRef + locaList[4]
        fixList = [fixChrom, fixStart, fixStop, fixRef, fixAlt]
    elif locaList[3] == "" and (locaList[4] != "" and locaList[4] != "-"):
        fixStart = str(int(locaList[1]) - 1)
        fixStop = locaList[2]
        fixRef = fixReference(fixChrom, fixStart, fixStop, reference)
        fixAlt = fixRef[0] + locaList[4]
        fixList = [fixChrom, fixStart, fixStop, fixRef, fixAlt]
    elif locaList[4] == "-":
        fixStart = str(int(locaList[1]) - 1)
        fixStop = locaList[2]
        fixRef = fixReference(fixChrom, fixStart, fixStop, reference)
        fixAlt = fixRef[0]
        fixList = [fixChrom, fixStart, fixStop, fixRef, fixAlt]
    else:
        fixList = locaList
    return fixList

# 将gDNA chr17:g.7577535C>G 转为annovar格式
o = open("gDNA.txt", "r", encoding="utf-8")
output = open("results.txt", "w", encoding="utf-8")
for line in o:
    line = line.replace("\n", "")
    l = fixFormat(line)
    output.write("\t".join(l) + "\n")
output.close()
o.close()


