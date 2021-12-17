# coding=utf-8
# pzw
# 20211216
# refTranscript
# MANE select > LRG > refSeq > Clinvar > HGNC

import gzip
from datetime import datetime
import sys

# MANE select
# https://ftp.ncbi.nlm.nih.gov/refseq/MANE/MANE_human/release_0.95/MANE.GRCh38.v0.95.summary.txt.gz
def maneDict(maneFile, refDict):
    mane = gzip.open(maneFile, "rb")
    for line in mane:
        line = line.decode("utf-8")
        if not line.startswith("#"):
            lines = line.split("\t")
            symbol = lines[3]
            transcript = lines[5]
            if not symbol in refDict.keys():
                refDict[symbol] = [transcript, "MANE"]
    mane.close()
    return refDict

# LRG
# ftp://ftp.ebi.ac.uk/pub/databases/lrgex/list_LRGs_transcripts_xrefs.txt
def lrgDict(lrgFile, refDict):
    lrg = open(lrgFile, "r", encoding="utf-8")
    for line in lrg:
        if not line.startswith("#"):
            lines = line.split("\t")
            symbol = lines[1]
            transcript = lines[4]
            if not symbol in refDict.keys():
                refDict[symbol] = [transcript, "LRG"]
            # LRG与MANE中倾向于更小的编号
            else:
                if int(transcript.split("_")[1].split(".")[0]) < int(refDict[symbol][0].split("_")[1].split(".")[0]):
                    refDict[symbol] = [transcript, "LRG"]
    lrg.close()
    return refDict

# refSeq
# https://ftp.ncbi.nlm.nih.gov/refseq/H_sapiens/RefSeqGene/LRG_RefSeqGene
def refseqDict(refseqFile, refDict):
    refseq = open(refseqFile, "r", encoding="utf-8")
    for line in refseq:
        if not line.startswith("#"):
            lines = line.replace("\n", "").split("\t")
            if "reference" in lines[6]:
                symbol = lines[2]
                transcript = lines[4]
                if not symbol in refDict.keys():
                    refDict[symbol] = [transcript, "refSeq"]
    refseq.close()
    return refDict

# clinvar
# https://ftp.ncbi.nlm.nih.gov/pub/clinvar/tab_delimited/variant_summary.txt.gz
def clinvarDict(clinvarFile, refDict):
    clinvar = gzip.open(clinvarFile, "rb")
    for line in clinvar:
        line = line.decode("utf-8")
        if not line.startswith("#"):
            lines = line.split("\t")
            symbol = lines[4]
            try:
                transcript = lines[3].split(":")[0].split("(")[0]
            except:
                transcript = "-"
            if "NM" in transcript:
                if not symbol in refDict.keys():
                    refDict[symbol] = [transcript, "Clinvar"]
                else:
                    if refDict[symbol][1] == "Clinvar":
                        if int(transcript.split("_")[1].split(".")[0]) < int(refDict[symbol][0].split("_")[1].split(".")[0]):
                            refDict[symbol] = [transcript, "Clinvar"]
    clinvar.close()
    return refDict

# HGNC
# https://ftp.ebi.ac.uk/pub/databases/genenames/hgnc/tsv/locus_types/gene_with_protein_product.txt
def hgncDict(hgncFile, refDict):
    hgnc = open(hgncFile, "r", encoding="utf-8")
    for line in hgnc:
        if not line.startswith("hgnc_id"):
            lines = line.split("\t")
            symbol = lines[1]
            transcript = lines[52]
            if transcript != "":
                transcript = transcript.replace('"', '').split("|")[1]
                if not symbol in refDict.keys():
                    refDict[symbol] = [transcript, "HGNC"]
            else:
                transcript = lines[23]
                if transcript != "":
                    if not symbol in refDict.keys():
                        refDict[symbol] = [transcript, "HGNC"]
    hgnc.close()
    return refDict

# APPEND
def main(database):
    refDict = {}
    maneDict(database + "/MANE.GRCh38.v0.95.summary.txt.gz", refDict)
    lrgDict(database + "/list_LRGs_transcripts_xrefs.txt", refDict)
    refseqDict(database + "/LRG_RefSeqGene", refDict)
    clinvarDict(database + "/variant_summary.txt.gz", refDict)
    hgncDict(database + "/gene_with_protein_product.txt", refDict)

    date = datetime.now().strftime("%Y%m%d")

    refTranscript = open("refTranscript_" + date + ".txt", "w", encoding="utf-8")
    refTranscript.write("#Gene\tTranscript\tSource\n")
    l = sorted(list(refDict.keys()))
    for k in l:
        refDict[k].insert(0, k)
        refTranscript.write("\t".join(refDict[k]) + "\n")
    refTranscript.close()

main(sys.argv[1])
