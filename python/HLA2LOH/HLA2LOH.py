# coding=utf-8
# pzw
# HLA2LOH
# 分析HLA Class I 区域 LOH

import os

# 提取normal与tumor的HLA区域reads
# 依赖samtools, sambamba
def extractHLAReads(bam, sampleID, buildver, threads, outputDir):
    # HLA 提取区域
    if buildver == "hg19":
        HLA = [
            "chr6:28477797-33448354",
            "chr6_apd_hap1:1-4622290",
            "chr6_cox_hap2:1-4795371",
            "chr6_dbb_hap3:1-4610396",
            "chr6_mann_hap4:1-4683263",
            "chr6_mcf_hap5:1-4833398",
            "chr6_qbl_hap6:1-4611984",
            "chr6_ssto_hap7:1-4928567"
        ]

    elif buildver == "hg38":
        HLA = [
            "chr6:28510120-33480577",
            "chr6_GL000250v2_alt:1066038-4433734",
            "chr6_GL000251v2_alt:1283988-4540572",
            "chr6_GL000252v2_alt:1063230-4372611",
            "chr6_GL000253v2_alt:1062914-4548533",
            "chr6_GL000254v2_alt:1062887-4416229",
            "chr6_GL000255v2_alt:1063190-4323464",
            "chr6_GL000256v2_alt:1106450-4577757"
        ]
    else:
        print("Cannot get buildver")
        exit()

    print("正在提取HLA区域reads")
    extractRegion = " ".join(HLA)
    if not os.path.exists(outputDir):
        os.makedirs(outputDir)
    if not os.path.exists(outputDir + "/bam_" + sampleID):
        os.makedirs(outputDir + "/bam_" + sampleID)
    if not os.path.exists(outputDir + "/reads"):
        os.makedirs(outputDir + "/reads")
    
    cmd = """
        samtools view {bam} \\
            {extractRegion} \\
            -b > {outputDir}/bam_{sampleID}/{sampleID}.HLA.bam
        samtools view {bam} -bh -f 12 -@ {threads} > {outputDir}/bam_{sampleID}/{sampleID}.unmapped.bam
        samtools merge {outputDir}/bam_{sampleID}/{sampleID}.merge.bam {outputDir}/bam_{sampleID}/{sampleID}.HLA.bam {outputDir}/bam_{sampleID}/{sampleID}.unmapped.bam -@ {threads}
        samtools sort {outputDir}/bam_{sampleID}/{sampleID}.merge.bam -n -@ {threads} -o {outputDir}/bam_{sampleID}/{sampleID}.sort.bam
        samtools fastq {outputDir}/bam_{sampleID}/{sampleID}.sort.bam \\
            -1 {outputDir}/reads/{sampleID}.HLA.R1.fastq \\
            -2 {outputDir}/reads/{sampleID}.HLA.R2.fastq -@ {threads} -s /dev/null
        rm -rf {outputDir}/bam_{sampleID}
    """.format(bam=bam, extractRegion=extractRegion, outputDir=outputDir, sampleID=sampleID, threads=threads)
    os.system(cmd)

# 获得HLA分型
# 依赖optitype
def HLATyping(bam, normalID, outputDir, threads, buildver):
    extractHLAReads(bam, normalID, buildver, threads, outputDir)
    print("使用OptiType进行HLA分型")
    tmpDir = outputDir + "/optitype"
    if not os.path.exists(tmpDir):
        os.makedirs(tmpDir)
    optitype = "/home/bioinfo/ubuntu/software/OptiType-1.3.5/OptiTypePipeline.py"
    cmd = """
        python {optitype} \\
            -i {outputDir}/reads/{normalID}.HLA.R1.fastq {outputDir}/reads/{normalID}.HLA.R2.fastq \\
            -d -o {tmpDir} -p {normalID} -v
        cp {tmpDir}/{normalID}_result.tsv {outputDir}/{normalID}.HLA.results.txt
    """.format(optitype=optitype, tmpDir=tmpDir, normalID=normalID, outputDir=outputDir)
    os.system(cmd)
    print("HLA分型完成")

# 获取HLA ref.fa
# 数据库来源polysolver
# https://github.com/jason-weirather/hla-polysolver/blob/master/data/abc_complete.fasta
# 依赖bwa、sambamba、samtools
def extractHLAFasta(outputDir, normalID):
    print("提取HLA分型对应序列形成参考基因组")
    tmpDir = outputDir + "/ref"
    if not os.path.exists(tmpDir):
        os.makedirs(tmpDir)

    # 初始化HLA数据库
    db_dict = {}
    db = open("/home/bioinfo/ubuntu/software/HLA2LOH/databases/abc_complete.fasta", "r")
    for d in db:
        if ">hla" in d:
            hla_type = d[1:12]
            db_dict[hla_type] = []
        else:
            db_dict[hla_type].append(d)
    db.close()

    # 获得HLA结果
    hlas = open(outputDir + "/optitype/" + normalID + "_result.tsv", "r")
    HLA_list = []
    for line in hlas:
        if not "Objective" in line:
            if not line == "\n":
                lines = line.split("\t")
                HLA_list.append(lines[1][0:7].replace("A*", "hla_a_").replace(":", "_"))
                HLA_list.append(lines[2][0:7].replace("A*", "hla_a_").replace(":", "_"))
                HLA_list.append(lines[3][0:7].replace("B*", "hla_b_").replace(":", "_"))
                HLA_list.append(lines[4][0:7].replace("B*", "hla_b_").replace(":", "_"))
                HLA_list.append(lines[5][0:7].replace("C*", "hla_c_").replace(":", "_"))
                HLA_list.append(lines[6][0:7].replace("C*", "hla_c_").replace(":", "_"))
    hlas.close()

    # 获得HLA reference
    outputFasta = open(tmpDir + "/HLA.fasta", "w")
    outputBed = open(tmpDir + "/HLA.bed", "w")
    for h in HLA_list:
        outputFasta.write(">" + h + "\n")
        outputFasta.write(db_dict[h][0])
        end = str(len(db_dict[h][0]))
        start = "1"
        chrom = h
        outputBed.write("\t".join([chrom, start, end]) + "\n")
    outputFasta.close()
    outputBed.close()
    cmd = """
        bwa index {tmpDir}/HLA.fasta
        samtools faidx {tmpDir}/HLA.fasta
    """.format(tmpDir=tmpDir)
    os.system(cmd)

# 重新比对到HLA reference
def ReMapping(sampleID, threads, outputDir):
    print(sampleID + " 重新比对")
    if not os.path.exists(outputDir + "/bam"):
        os.makedirs(outputDir + "/bam")
    cmd = """
        bwa mem -t {threads} -M \\
            -R "@RG\\tSM:{sampleID}\\tID:{sampleID}\\tPU:LOH\\tPL:illumina" \\
            {outputDir}/ref/HLA.fasta \\
            {outputDir}/reads/{sampleID}.HLA.R1.fastq {outputDir}/reads/{sampleID}.HLA.R1.fastq \\
            | sambamba view -S -f bam /dev/stdin -t {threads} > {outputDir}/bam/{sampleID}.bam
        sambamba sort {outputDir}/bam/{sampleID}.bam -t {threads} -o {outputDir}/bam/{sampleID}.sort.bam
        rm {outputDir}/bam/{sampleID}.bam
        sambamba markdup -t {threads} {outputDir}/bam/{sampleID}.sort.bam {outputDir}/bam/{sampleID}.marked.bam
        rm {outputDir}/bam/{sampleID}.sort.bam*
        samtools view -bh -F 1280 {outputDir}/bam/{sampleID}.marked.bam -o {outputDir}/bam/{sampleID}.aln.bam
        samtools index {outputDir}/bam/{sampleID}.aln.bam
        rm {outputDir}/bam/{sampleID}.marked.bam*
    """.format(threads=threads, sampleID=sampleID, outputDir=outputDir)
    os.system(cmd)

# 计算LOH
# 依赖PyLOH







normalbam = "/home/bioinfo/ubuntu/project/LOH_20210423/bam/W22805.bam"
normalID = "W22805"
tumorbam = "/home/bioinfo/ubuntu/project/LOH_20210423/bam/T22805.bam"
tumorID = "T22085"
buildver = "hg19"
threads = "8"
outputDir = "/home/bioinfo/ubuntu/project/LOH_20210423/output"


# extractHLAReads(normalbam, normalID, buildver, threads, outputDir)
# HLATyping(normalbam, normalID, outputDir, threads, buildver)
# extractHLAFasta(outputDir, normalID)
# extractHLAReads(tumorbam, tumorID, buildver, threads, outputDir)
# ReMapping(normalID, threads, outputDir)
# ReMapping(tumorID, threads, outputDir)

