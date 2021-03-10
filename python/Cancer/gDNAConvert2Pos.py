# coding=utf-8
# pzw
# 20210310

import re
import sys


def gDNA2Pos(gdna):
    if (gdna == "") or (gdna == "-"):
        chrom = ""
        start = ""
        end = ""
        ref = ""
        alt = ""

    else:
        gdnas = gdna.split(":g.")
        chrom = gdnas[0]
        seqs = gdnas[1]
        if "_" in seqs:
            if "(" in seqs:
                # chr9:g.(130872876_130872877)
                seqss = seqs.replace("(", "").replace(")", "").split("_")
                start = seqss[0]
                end = seqss[1]
                ref = ""
                alt = ""

            else:
                if ("del" in seqs) and ("ins" in seqs):
                    # chr9:g.130872141_130872142delGAinsTG
                    seqss = seqs.split("_")
                    start = seqss[0]
                    seqsss = seqss[1].split("del")
                    end = seqsss[0]
                    seqssss = seqsss[1].split("ins")
                    ref = seqssss[0]
                    alt = seqssss[1]
                else:
                    if "del" in seqs:
                        # chr5:g.112838448_112838449delGC
                        seqss = seqs.split("_")
                        start = seqss[0]
                        seqsss = seqss[1].split("del")
                        end = seqsss[0]
                        ref = seqsss[1]
                        if ref.isdigit():
                            ref = ref + "bp"
                        alt = "-"
                    elif "ins" in seqs:
                        # chr5:g.112839932_112839933insGACATCGCA
                        seqss = seqs.split("_")
                        start = seqss[0]
                        seqsss = seqss[1].split("ins")
                        end = seqsss[0]
                        ref = "-"
                        alt = seqsss[1]
                    else:
                        seqss = seqs.split("_")
                        start = seqss[0]
                        end = seqss[1]
                        ref = ""
                        alt = ""                 

        else:
            if "dup" in seqs:
                # chr3:g.128481912dupT
                seqss = seqs.split("dup")
                start = seqss[0]
                end = seqss[0]
                ref = "-"
                alt = seqss[1]
            elif "del" in seqs:
                # chr3:g.128486343delG
                seqss = seqs.split("del")
                start = seqss[0]
                end = seqss[0]
                ref = seqss[1]
                alt = "-"
            elif ">" in seqs:
                # chr11:g.108272782G>T
                seqss = seqs.split(">")
                alt = seqss[1]
                splitNum = re.findall(r"[A-Za-z]+|\d+", seqss[0])
                start = splitNum[0]
                ref = splitNum[1]
                end = str(int(start) + len(ref) - 1)
            else:
                start = ""
                end = ""
                ref = ""
                alt = ""

    return [chrom, start, end, ref, alt]


a = open("gDNA.txt", "r")
b = open("results.txt", "w")
for line in a:
    lines = line.split("\n")[0]
    l = gDNA2Pos(lines)
    b.write("\t".join(l) + "\n")
b.close()
a.close()
