# coding=utf-8
# pzw
# NIPT analysis
# 20211101
# 依赖于bwa, samtools, sambamba, bedtools, R

import os
import sys
from numpy import sqrt
import pandas as pd
import argparse
import shutil

# 获得程序运行路径
def getAbsPath():
    now = os.path.abspath(os.path.dirname(sys.argv[0]))
    return now

# 获得Contig大小，只需要常染色体chr1-chr22
# reference需要预先建立索引
# samtools faidx reference.fa
def GetFastaBinsGCBias(reference, tmpDir, output):
    if not os.path.exists(tmpDir):
        os.makedirs(tmpDir)
    fai = open(reference + ".fai", "r")
    contig = open(tmpDir + "/contig.size", "w")
    contigList = ["chr" + str(x) for x in range(1, 23)]
    for line in fai:
        lines = line.split("\t")
        chrom = lines[0]
        length = lines[1]
        if chrom in contigList:
            contig.write(chrom + "\t" + length + "\n")
    contig.close()
    fai.close()
    cmd = """
        bedtools makewindows -g {tmpDir}/contig.size -w 300000 > {tmpDir}/contig.300kb.bed
        bedtools nuc -fi {reference} -bed {tmpDir}/contig.300kb.bed | cut -f 1-3,5 > {output}
    """.format(reference=reference, tmpDir=tmpDir, output=output)
    os.system(cmd)

# 获得每个bam文件的每个bin的counts数
# bam文件来源
# bwa mem
# samtools sort
# sambamba markdup -r
def GetBamCounts(bamPath, gcBed, output):
    if not os.path.exists(output):
        os.makedirs(output)
    for bam in os.listdir(bamPath):
        if bam.endswith(".bam"):
            name = bam.replace(".bam", ".counts")
            bamFile = bamPath + "/" + bam
            cmd = """
                bedtools coverage -a {gcBed} -b {bamFile} | cut -f 1-5 \\
                    | sed '1i chr\\tstart\\tend\\tGC\\tRC' > {output}/{name}
            """.format(gcBed=gcBed, bamFile=bamFile, output=output, name=name)
            os.system(cmd)

# 适用GC比例校正read counts
def LoessCorrectCountsWithGC(countPath, output):
    if not os.path.exists(output):
        os.makedirs(output)
    script = getAbsPath() + "/gcContent.R"
    for count in os.listdir(countPath):
        gcCount = count.replace(".count", ".cor.count")
        cmd = """
            Rscript {script} {countPath}/{count} {output}/{gcCount}
        """.format(script=script, countPath=countPath, count=count, output=output, gcCount=gcCount)
        os.system(cmd)

# 形成Counts矩阵
def CountsMatrix(corPath, matrixFile):
    countDict = {}
    for txt in os.listdir(corPath):
        name = txt.replace(".cor.count", "")
        countDict[name] = {}
        f = open(corPath + "/" + txt, "r")
        for line in f:
            if not line.startswith("chr\tstart"):
                lines = line.replace("\n", "").split("\t")
                chromosome = lines[0]
                readCounts = lines[4]
                try:
                    countDict[name][chromosome] += float(readCounts)
                except:
                    countDict[name][chromosome] = 0
                    countDict[name][chromosome] += float(readCounts)
        f.close()
    matrix = open(matrixFile, "w")
    kList = list(countDict.keys())
    matrix.write("#SampleName\t" + "\t".join(kList) + "\n")
    for i in range(22):
        chr = "chr" + str(i+1)
        countList = []
        for k in kList:
            countList.append("%.2f" % countDict[k][chr])
        matrix.write(chr + "\t" + "\t".join(countList) + "\n")
    matrix.close()

# 计算referenceSD 以及 referenceMean
def URMatrix(countsMatrix, ur="reference.score"):
    df = pd.read_table(countsMatrix, sep="\t", header=0, index_col=0)
    colSumDict = {}
    for col in df.columns.values.tolist():
        colSumDict[col] = df[col].sum()
    df_UR = df
    for row in df_UR.index:
        for col in colSumDict.keys():
            df_UR.loc[row, col] = df_UR.loc[row, col] / colSumDict[col]
    rowSumDict = {}
    for row in df_UR.index:
        rowSumDict[row] = 0
        for col in df_UR.columns.values.tolist():
            rowSumDict[row] += df_UR.loc[row, col]
    meanDict = {}
    sampleAmount = df_UR.shape[1]
    for chrom in rowSumDict.keys():
        meanDict[chrom] = rowSumDict[chrom] / sampleAmount
    sdDict = {}
    for row in df_UR.index:
        sqSum = 0
        for col in df_UR.columns.values.tolist():
            sqSum += (df_UR.loc[row, col] - meanDict[row]) ** 2
        sdDict[row] = sqrt(sqSum / sampleAmount)
    
    ref = open(ur, "w")
    ref.write("#Chrom\tMean\tSD\n")
    for k in meanDict.keys():
        meanChr = meanDict[k]
        sdChr = sdDict[k]
        ref.write(k + "\t" + str(meanChr) + "\t" + str(sdChr) + "\n")
    ref.close()

def main(reference, bamPath, countsFile, refFile):
    tmpdir = getAbsPath() + "/tmp"
    GetFastaBinsGCBias(reference, tmpdir, tmpdir + "/gc.bed")
    GetBamCounts(bamPath, tmpdir + "/gc.bed", tmpdir + "/referenceCounts")
    LoessCorrectCountsWithGC(tmpdir + "/referenceCounts", tmpdir + "/referenceCorCounts")
    CountsMatrix(tmpdir + "/referenceCorCounts", countsFile)
    URMatrix(countsFile, refFile)
    shutil.rmtree(tmpdir)
    print("成功生成参考文件：", refFile)

if __name__ == "__main__":
    parser = argparse.ArgumentParser(
        description="make NIPT reference @pzweuj",
        prog="reference.py",
        usage="python3 reference.py [-h] -f <ref.fa> -b <bam files dir path> -o <ref file> -c <counts output file>",
        formatter_class=argparse.RawTextHelpFormatter
    )
    parser.add_argument("-v", "--version", action="version",
        version="Version 0.1 20211104")
    parser.add_argument("-f", "--fasta", type=str,
        help="reference.fasta(contain chr1-22)")
    parser.add_argument("-b", "--bamPath", type=str,
        help="reference bam files directory")
    parser.add_argument("-c", "--counts", type=str,
        help="read counts output")
    parser.add_argument("-o", "--output", type=str,
        help="z score ref file")

    if len(sys.argv[1:]) == 0:
        parser.print_help()
        parser.exit()
    args = parser.parse_args()
    main(reference=args.fasta, bamPath=args.bamPath, countsFile=args.counts, refFile=args.output)
