# coding=utf-8
# 20210204
# pzw
# 体细胞变异分析流程
# 参考基因组：hg19

import os
import time
import json


## 时间获取
def getTime():
	localTime = time.asctime(time.localtime(time.time()))
	return localTime

## 参考基因组
### hg19
### 使用hg19以及相关的数据库
def ReferenceVersion(Ref="hg19"):
	if Ref == "hg19":
		ReferenceDict = {
			"reference": "hg19.fasta",
			"indelReference1": "Mills_and_1000G_gold_standard.indels.hg19.sites.vcf",
			"indelReference2": "1000G_phase1.indels.hg19.sites.vcf",
			"gnomad": "af-only-gnomad.raw.sites.hg19.vcf.gz",
			"humandb": "/humandb_hg19/",
			"mutect2_pon": "pon_dir",
			"pisces_reference": "/hg19"
		}
	else:
		ReferenceDict = {}
		print("默认基因组hg19，暂未设定其他参考！")
		exit()
	return ReferenceDict

## 样本信息预处理
### 使用json作为参数输入文件
### 包括样本信息的获取，相关输出文件夹的建立等
def jsonFileReader(jsonFile):
	data = json.dump(jsonFile)
	return data

## 质控流程
### fastp质控
### 如有确认的UMI，调整模块
def fastpFilter(sampleID, rawdataDir, outputDir, threads, UMI):
	print(sampleID + " Fastp过滤开始", getTime())
	cmd = """
		fastp -i {rawdataDir}/{sampleID}.R1.fastq.gz \\
			-I {rawdataDir}/{sampleID}.R2.fastq.gz \\
			-o {outputDir}/{sampleID}.clean.R1.fastq.gz \\
			-O {outputDir}/{sampleID}.clean.R2.fastq.gz \\
			-w {threads} -j {outputDir}/{sampleID}.json \\
			-h {outputDir}/{sampleID}.html --umi_len {UMI}
	""".format(rawdataDir=rawdataDir, sampleID=sampleID, outputDir=outputDir, threads=threads, UMI=UMI)
	os.system(cmd)
	print(sampleID + " Fastp过滤结束", getTime())


## 比对流程
### bwa比对
def bwa_mapping(sampleID, cleandataDir, outputDir, threads):
	print(sampleID + " BWA比对开始", getTime())
	reference = ReferenceVersion("hg19")["reference"]
	cmd = """
		bwa mem -t {threads} {reference} \\
			{cleandataDir}/{sampleID}.clean.R1.fastq.gz \\
			{cleandataDir}/{sampleID}.clean.R2.fastq.gz \\
			| samtools view -bSh - \\
			| samtools sort -@ {threads} - -o {outputDir}/{sampleID}.bam
	""".format(threads=threads, cleandataDir=cleandataDir, outputDir=outputDir, reference=reference)
	os.system(cmd)
	print(sampleID + " BWA比对结束", getTime())


## 重校正流程
### 使用GATK4对indel区域进行重校正
### 最终保留原始比对的bam和校正后的bam
def gatk4_recalibrator(sampleID, tempBamDir, finalBamDir):
	print(sampleID + " GATK4重校正开始", getTime())
	reference = ReferenceVersion("hg19")["reference"]
	indelReference1 = ReferenceVersion("hg19")["indelReference1"]
	indelReference2 = ReferenceVersion("hg19")["indelReference2"]
	cmd = """
		gatk MarkDuplicates -I {finalBamDir}/{sampleID}.bam \\
			-O {tempBamDir}/{sampleID}.mdups.bam
		gatk AddOrReplaceReadGroups -I {tempBamDir}/{sampleID}.mdups.bam \\
			-O {tempBamDir}/{sampleID}.header.bam \\
			-LB somatic -PL illumina -PU somatic -SM {sampleID}
		gatk BuildBamIndex -I {tempBamDir}/{sampleID}.header.bam
		gatk BaseRecalibrator -I {tempBamDir}/{sampleID}.header.bam \\
			-R {reference} \\
			--known-sites {indelReference1} \\
			--known-sites {indelReference2} \\
			-O {tempBamDir}/{sampleID}.recal.table
		gatk ApplyBQSR -R {reference} \\
			-I {tempBamDir}/{sampleID}.header.bam \\
			--bqsr {tempBamDir}/{sampleID}.recal.table \\
			-O {finalBamDir}/{sampleID}.BQSR.bam
	""".format(tempBamDir=tempBamDir, sampleID=sampleID, reference=reference, indelReference1=indelReference1, indelReference2=indelReference2)
	os.system(cmd)
	print(sampleID + " GATK4重校正结束", getTime())


## SNV/indel变异检测
### 文章DOI: 10.1038/srep43169测试在VAF<=0.10时，表现最佳软件为VarDict，在VAF>0.10时，GATK表现较好。
### 使用GATK4 HaplotypeCaller进行检测
def gatk4_haplotypecaller_caller(sampleID, finalBamDir, VcfDir, threads):
	print(sampleID, " HaplotypeCaller检测开始", getTime())
	reference = ReferenceVersion("hg19")["reference"]
	cmd = """
		gatk HaplotypeCaller \\
			-R {reference} \\
			-I {finalBamDir}/{sampleID}.BQSR.bam \\
			--native-pair-hmm-threads {threads} \\
			-O {VcfDir}/{sampleID}.haplotypeCaller.vcf
	""".format(reference=reference, finalBamDir=finalBamDir, sampleID=sampleID, threads=threads)
	os.system(cmd)
	print(sampleID, " HaplotypeCaller检测结束", getTime())

### 使用GATK4 Mutect2进行检测
#### Mutect2单样本模式要求使用多个正常样本构建基线
#### 详见https://gatk.broadinstitute.org/hc/en-us/articles/360046788432-Mutect2
"""
gatk Mutect2 -R reference.fasta -I normal.bam -max-mnp-distance 0 \
	-O normal.vcf.gz
gatk GenomicsDBImport -R referencefasta -L intervals.interval_list \
	--genomicsdb-workspace-path pon-db \
	-V normal1.vcf.gz \
	-V normal2.vcf.gz \
	-V normal3.vcf.gz
gatk CreateSomaticPanelOfNormals -R reference.fasta -V pon_db \
	-O pon.vcf.gz
"""
def gatk4_mutect2_caller(sampleID, finalBamDir, VcfDir, threads):
	print(sampleID, " Mutect2检测开始", getTime())
	reference = ReferenceVersion("hg19")["reference"]
	germline_resource = ReferenceVersion("hg19")["gnomad"]
	pon = ReferenceVersion("hg19")["mutect2_pon"]
	cmd = """
		gatk Mutect2 -R {reference} \\
			-I {finalBamDir}/{sampleID}.BQSR.bam \\
			--germline-resource {germline_resource} \\
			--panel-of-normals {pon} \\
			--native-pair-hmm-threads {threads}
			-O {VcfDir}/{sampleID}.mutect2.vcf
	""".format(reference=reference, finalBamDir=finalBamDir, sampleID=sampleID, germline_resource=germline_resource, pon=pon, threads=threads, VcfDir=VcfDir)
	os.system(cmd)
	print(sampleID, " Mutect2检测结束", getTime())

### 使用freebayes进行检测
def freebayes_caller(sampleID, finalBamDir, VcfDir, threads):
	print(sampleID, " freebayes检测开始", getTime())
	reference = ReferenceVersion("hg19")["reference"]
	cmd = """
		freebayes -f {reference} {finalBamDir}/{sampleID}.BQSR.bam \\
			> {VcfDir}/{sampleID}.freebayes.vcf
	""".format(reference=reference, finalBamDir=finalBamDir, sampleID=sampleID, VcfDir=VcfDir)
	os.system(cmd)
	print(sampleID, " freebayes检测结束", getTime())

### 使用Pisces进行检测
### 使用Pisces需先对参考基因组创建拥有Pisces的索引
"""
dotnet CreateGenomeSizeFile.dll \
	-s "Homo_Sapiens(UCSC hg19)" \
	-g /mnt/d/ubuntu/database/hg19/ \
	-out /mnt/d/ubuntu/database/hg19_pisces

# 以上目前测试有bug
# 详见https://github.com/Illumina/Pisces/issues/23
"""
def pisces_caller(sampleID, finalBamDir, VcfDir, threads)
	print(sampleID, " pisces检测开始", getTime())
	reference = ReferenceVersion("hg19")["pisces_reference"]
	cmd = """
		dotnet Pisces.dll -bam {finalBamDir}/{sampleID}.BQSR.bam \\
			-g {reference} --callmnvs false \\
			-gVCF false -o {VcfDir} -t {threads}
	""".format(finalBamDir=finalBamDir, sampleID=sampleID, reference=reference, VcfDir=VcfDir, threads=threads)
	os.system(cmd)
	print(sampleID, " pisces检测结束", getTime())


### 使用VarDict进行检测
### 设定检测阈值为0.002，ctDNA的检测限
def vardict_caller(sampleID, finalBamDir, VcfDir):
	print(sampleID, " vardict检测开始", getTime())
	reference = ReferenceVersion("hg19")["reference"]
	cmd = """
		perl vardict.pl -G {reference} -f 0.002 \\
			-N {sampleID} -b {finalBamDir}/{sampleID}.BQSR.bam \\
			| teststrandbias.R \\
			| var2vcf_valid.pl -N {VcfDir}/{sampleID} \\
			-E -f 0.002
	""".format(sampleID=sampleID, finalBamDir=finalBamDir, reference=reference, VcfDir=VcfDir)
	os.system(cmd)
	print(sampleID, " vardict检测结束", getTime())

## 过滤
### 使用bedtools根据bed文件进行过滤
def bedtools_filter(sampleID, VcfDir, bedFile):
	print(sampleID, "bedtools过滤开始", getTime())
	cmd = """
		bedtools intersect -a {VcfDir}/{sampleID}.vcf -b {bedFile} \
			> {VcfDir}/{sampleID}.target.vcf
	""".format(VcfDir=VcfDir, sampleID=sampleID, bedFile=bedFile)
	os.system(cmd)
	print(sampleID, "bedtools过滤结束", getTime())

## 注释
### 使用annovar进行注释
### 测试中annovar数据库未实装，实际使用数据库需调整
def annovar_anno(sampleID, VcfDir, AnnotationDir, threads):
	print(sampleID, " annovar注释开始", getTime())
	humandb = ReferenceVersion("hg19")["humandb_hg19"]
	cmd = """
		convert2annovar.pl -format vcf4 {VcfDir}/{sampleID}.vcf > {AnnotationDir}/{sampleID}.avinput
		table_annovar.pl {AnnotationDir}/{sampleID}.avinput \\
			{humandb} -buildver hg19 \\
			-out {AnnotationDir}/{sampleID} -remove \\
			-protocol refGene,cytoBand,avsnp150,gnomad211,clinvar,cosmic,dbnsfp35a \\
			-operation g,r,f,f,f,f,f \\
			-nastring . -thread {threads} -otherinfo
	""".format(VcfDir=VcfDir, sampleID=sampleID, AnnotationDir=AnnotationDir, humandb=humandb, threads=threads)
	os.system(cmd)
	print(sampleID, " annovar注释结束", getTime())

### 使用SnpEff进行注释
### SnpEff注释后生产vcf文件，可继续用于annovar的注释
def snpeff_anno(sampleID, VcfDir, AnnotationDir):
	print(sampleID + " snpEff注释开始", getTime())
	cmd = """
		java -jar snpEff.jar GRCh37 {VcfDir}/{sampleID}.vcf > {AnnotationDir}/{sampleID}.snpeff.vcf
	""".format(VcfDir=VcfDir, sampleID=sampleID, AnnotationDir=AnnotationDir)
	os.system(cmd)
	print(sampleID + " snpEff注释结束", getTime())

## MSI检测
### 使用MSIsensor2进行检测
def msisensor2_msi(bamFile, MSIDir):
	print(sampleID + " msisensor2检测开始", getTime())
	cmd = """
		msisensor2 msi -M models_hg19 -t {bamFile} -o {MSIDir}
	""".format(bamFile=bamFile, MSIDir=MSIDir)
	os.system(cmd)
	print(sampleID + " msisensor2检测结束", getTime())

## CNV检测
### Cnvkit检测
"""
python cnvkit.py reference -o FlatReference.cnn \
	-f reference.fa -t my_targets.bed -a my_antitargets.bed
python cnvkit.py fix sample.targetcoverage.cnn \
	sample.antitargetcoverage.cnn \
	reference.cnn \
	-o sample.cnr
python cnvkit.py segment \
	sample.cnr \
	-o sample.cns
python cnvkit.py call \
	sample.cns \
	-o sample.call.cns
"""


## SV检测
### lumpy检测
### lumpy检测需求使用samtools获得两个中间bam文件
def lumpy_sv(sampleID, finalBamDir, SvDir, threads):
	print(sampleID + " lumpy检测开始", getTime())
	cmd = """
		samtools view {finalBamDir}/{sampleID}.BQSR.bam \\
			-bh -F 1294 \\
			| samtools sort - -o {SvDir}/{sampleID}.discordants.bam
		samtools view -h {finalBamDir}/{sampleID}.BQSR.bam \\
			| extractSplitReads_BwaMem \\
			-i stdin \\
			| samtools view -bSh - \\
			| samtools sort -@ {threads} - -o {SvDir}/{sampleID}.splitters.bam
		lumpyexpress -B {finalBamDir}/{sampleID}.BQSR.bam \\
			-D {SvDir}/{sampleID}.discordants.bam \\
			-S {SvDir}/{sampleID}.splitters.bam \\
			-o {SvDir}/{sampleID}.lumpy.vcf
		svtyper-sso -i {SvDir}/{sampleID}.lumpy.vcf \\
			-B {finalBamDir}/{sampleID}.BQSR.bam \\
			--cores {threads} \\
			-o {SvDir}/{sampleID}.svtyper.vcf
	""".format(finalBamDir=finalBamDir, sampleID=sampleID, SvDir=SvDir, threads=threads)
	os.system(cmd)
	print(sampleID + " lumpy结束", getTime())

### FACTERA检测
### FACTERA许可证声明禁止盈利使用
def factera_sv(sampleID, finalBamDir, SvDir):
	print(sampleID + " factera检测开始", getTime())
	reference = "hg19.2bit"
	cmd = """
		perl factera.pl {finalBamDir}/{sampleID}.BQSR.bam \\
			{reference} -o {SvDir}
	""".format(finalBamDir=finalBamDir, sampleID=sampleID, SvDir=SvDir, reference=reference)
	os.system(cmd)
	print(sampleID + " factera检测结束", getTime())