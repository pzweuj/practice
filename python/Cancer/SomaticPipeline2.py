# coding=utf-8
# 20210207
# pzw
# 体细胞变异分析流程
# 参考基因组：hg19

import os

"""
#软件目录 地址 开源许可
fastp https://github.com/OpenGene/fastp MIT
fgbio https://github.com/fulcrumgenomics/fgbio MIT
bwa https://github.com/lh3/bwa GPLv3
samtools https://github.com/samtools/samtools MIT
GATK4 https://github.com/broadinstitute/gatk BSD-3
freebayes https://github.com/freebayes/freebayes MIT
annovar https://annovar.openbioinformatics.org/en/latest/user-guide/download/ non-profit-use-only
cnvkit https://cnvkit.readthedocs.io/en/stable/pipeline.html Apache
lumpy https://github.com/arq5x/lumpy-sv MIT
svtyper https://github.com/hall-lab/svtyper MIT
msisensor2 https://github.com/niu-lab/msisensor2 GPLv3
fgbio http://fulcrumgenomics.github.io/fgbio/tools/latest/ MIT
"""


"""
#数据库 来源
ucsc.hg19.fasta gsapubftp-anonymous@ftp.broadinstitute.org:/bundle/hg19
Mills_and_1000G_gold_standard.indels.hg19.sites.vcf.gz gsapubftp-anonymous@ftp.broadinstitute.org:/bundle/hg19
dbsnp_138.hg19.vcf.gz gsapubftp-anonymous@ftp.broadinstitute.org:/bundle/hg19
1000G_phase1.snps.high_confidence.hg19.sites.vcf.gz gsapubftp-anonymous@ftp.broadinstitute.org:/bundle/hg19
1000G_phase1.indels.hg19.sites.vcf.gz gsapubftp-anonymous@ftp.broadinstitute.org:/bundle/hg19
af-only-gnomad.raw.sites.hg19.vcf.gz
small_exac_common_3_b37.vcf
illumina_pt2.interval_list
reFlat.txt http://hgdownload.soe.ucsc.edu/goldenPath/hg19/database/refFlat.txt.gz
wgEncodeDukeMapabilityRegionsExcludable.bed ftp://hgdownload.soe.ucsc.edu/goldenPath/hg19/encodeDCC/wgEncodeMapability/wgEncodeDukeMapabilityRegionsExcludable.bed.gz
"""

## 创建文件夹
def mkdir(folder_path):
	folder = os.path.exists(folder_path)
	if not folder:
		os.makedirs(folder_path)
		print("create new folder " + folder_path + " OK")
	else:
		print("folder " + folder_path + " exists")

## 质控流程
### 使用fastp
def fastp_qc(sampleID, rawdataDir, resultsDir, threads):
	mkdir(resultsDir)
	mkdir(resultsDir + "/QC")
	mkdir(resultsDir + "/cleandata")
	cmd = """
		fastp -i {rawdataDir}/{sampleID}_R1.fastq.gz \\
			-I {rawdataDir}/{sampleID}_R2.fastq.gz \\
			-o {resultsDir}/cleandata/{sampleID}.clean_R1.fastq.gz \\
			-O {resultsDir}/cleandata/{sampleID}.clean_R2.fastq.gz \\
			-j {resultsDir}/QC/{sampleID}.json \\
			-h {resultsDir}/QC/{sampleID}.html \\
			-w {threads}
	""".format(rawdataDir=rawdataDir, sampleID=sampleID, resultsDir=resultsDir, threads=threads)
	os.system(cmd)

## UMI处理
### 使用fgbio，待定


## 比对
### 使用bwa以及samtools
def bwa_mapping(sampleID, resultsDir, threads):
	mkdir(resultsDir + "/bam")
	reference = "ucsc.hg19.fasta"
	cmd = """
		bwa mem -t {threads} \\
			-M \\
			-R "@RG\\tID:{sampleID}\\tLB:{sampleID}\\tPL:illumina\\t:PU:Hiseq\\tSM:{sampleID}" \\
			{reference} \\
			{resultsDir}/cleandata/{sampleID}.clean_R1.fastq.gz \\
			{resultsDir}/cleandata/{sampleID}.clean_R2.fastq.gz \\
			| samtools view -bSh --threads {threads} - > {resultsDir}/bam/{sampleID}.bam
			samtools index {resultsDir}/bam/{sampleID}.bam
	""".format(threads=threads, sampleID=sampleID, reference=reference, resultsDir=resultsDir)
	os.system(cmd)

## 标记重复
### 使用GATK4
def gatk4_markdups(sampleID, resultsDir):
	mkdir(resultsDir + "/tempFile")
	cmd = """
		gatk MarkDuplicates \\
			-I {resultsDir}/bam/{sampleID}.bam \\
			-O {resultsDir}/tempFile/{sampleID}_marked.bam \\
			-M {resultsDir}/tempFile/{sampleID}_marked_metrics.txt
	""".format(resultsDir=resultsDir, sampleID=sampleID)
	os.system(cmd)

## 校正
### 使用GATK4
def gatk4_recal(sampleID, resultsDir):
	mkdir(resultsDir + "/tempFile")
	reference = "ucsc.hg19.fasta"
	dbsnp = "dbsnp_138.hg19.vcf.gz"
	mills = "Mills_and_1000G_gold_standard.indels.hg19.sites.vcf.gz"
	G1000 = "1000G_phase1.indels.hg19.sites.vcf.gz"
	cmd = """
		gatk BaseRecalibrator \\
			--known-sites {dbsnp} \\
			--known-sites {mills} \\
			--known-sites {G1000} \\
			-R {reference} \\
			-I {resultsDir}/tempFile/{sampleID}_marked.bam \\
			-O {resultsDir}/tempFile/{sampleID}_recal.table

		gatk ApplyBQSR \\
			-R {reference} \\
			--bqsr-recal-file {resultsDir}/tempFile/{sampleID}_recal.table \\
			-I {resultsDir}/tempFile/{sampleID}_marked.bam \\
			-O {resultsDir}/bam/{sampleID}_bqsr.bam
	""".format(dbsnp=dbsnp, mills=mills, G1000=G1000, reference=reference, resultsDir=resultsDir, sampleID=sampleID)
	os.system(cmd)


## 变异检测
### 使用GATK4
### GATK4需求先使用正常样本建立Panel of Normal
### 参考https://gatk.broadinstitute.org/hc/en-us/articles/360035531132--How-to-Call-somatic-mutations-using-GATK4-Mutect2
#### 其中normal bam可按tumor bam同样的流程生成
def gatk4_pon1(normalBam, normalVcf):
	reference = "ucsc.hg19.fasta"
	cmd = """
		gatk Mutect2 -R {reference} \\
			-I {normalBam} --max-mnp-distance 0 \\
			-O {normalVcf}
	""".format(reference=reference, normalBam=normalBam, normalVcf=normalVcf)
	os.system(cmd)

#### 把每个normal vcf使用germline数据库校对并存入pon数据中
def gatk4_pon2(normalVcf, pon_db):
	reference = "ucsc.hg19.fasta"
	cmd = """
		gatk GenomicsDBImport -R {reference} \\
			--genomicsdb-workspace-path {pon_db} \\
			-V {normalVcf}
	""".format(reference=reference, normalVcf=normalVcf, pon_db=pon_db)
	os.system(cmd)

#### 从pon数据库生成pon参考集
def gatk4_pon3(pon_db, pon):
	reference = "ucsc.hg19.fasta"
	gnomad = "af-only-gnomad.raw.sites.hg19.vcf.gz"
	cmd = """
		gatk CreateSomaticPanelOfNormals -R {reference} \\
			--germline-resource {gnomad} \\
			-V {pon_db} \\
			-O {pon}
	""".format(reference=reference, gnomad=gnomad, pon_db=pon_db, pon=pon)
	os.system(cmd)

#### 污染估算
#### 参考https://docs.gdc.cancer.gov/Data/Bioinformatics_Pipelines/DNA_Seq_Variant_Calling_Pipeline/#tumor-only-variant-calling-workflow
def gatk4_calCon():
	small_exac = "small_exac_common_3_b37.vcf"
	interList = "illumina_pt2.interval_list"
	reference = "ucsc.hg19.fasta"
	cmd = """
		gatk GetPileupSummaries \\
			-I {resultsDir}/bam/{sampleID}_bqsr.bam \\
			-O {resultsDir}/tempFile/{sampleID}_pileups.table \\
			-V {small_exac} \\
			-L {interList} \\
			-R {reference}

		gatk CalculateContamination \\
			-I {resultsDir}/tempFile/{sampleID}_pileups.table \\
			-O {resultsDir}/tempFile/{sampleID}_contamination.table
	""".format(resultsDir=resultsDir, sampleID=sampleID, small_exac=small_exac, interList=interList, reference=reference)
	os.system(cmd)

#### Mutect2变异检测
def gatk4_m2(sampleID, resultsDir, threads, pon, threads):
	mkdir(resultsDir + "/vcf")
	reference = "ucsc.hg19.fasta"
	gnomad = "af-only-gnomad.raw.sites.hg19.vcf.gz"
	cmd = """
		gatk Mutect2 \\
			-R {reference} \\
			-I {resultsDir}/bam/{sampleID}_bqsr.bam \\
			-O {resultsDir}/vcf/{sampleID}.m2.vcf \\
			-tumor {sampleID} \\
			--af-of-alleles-not-in-resource 0.0000025 \\
			--germline-resource {gnomad} \\
			-pon {pon} \\
			--native-pair-hmm-threads {threads}
	""".format(reference=reference, resultsDir=resultsDir, sampleID=sampleID, gnomad=gnomad, pon=pon, threads=threads)
	os.system(cmd)

#### Mutect2结果过滤
def gatk4_filter(resultsDir, sampleID):
	reference_dict = "ucsc.hg19.dict"
	cmd = """
		gatk SortVcf \\
			-I {resultsDir}/vcf/{sampleID}.m2.vcf \\
			-SD {reference_dict} \\
			-O {resultsDir}/tempFile/{sampleID}.m2.sorted.vcf

		gatk FilterMutectCalls \\
			-V {resultsDir}/tempFile/{sampleID}.m2.sorted.vcf \\
			-O {resultsDir}/vcf/{sampleID}.m2.contFiltered.vcf \\
			--contamination-table {resultsDir}/tempFile/{sampleID}_contamination.table
	""".format(resultsDir=resultsDir, sampleID=sampleID, reference_dict=reference_dict)
	os.system(cmd)

### 使用freebayes
def freebayes(resultsDir, sampleID):
	mkdir(resultsDir + "/vcf")
	reference = "ucsc.hg19.fasta"
	cmd = """
		freebayes -f {reference} \\
			{resultsDir}/bam/{sampleID}_bqsr.bam \\
			> {resultsDir}/vcf/{sampleID}.freebayes.vcf
	""".format(reference=reference, resultsDir=resultsDir, sampleID=sampleID)
	os.system(cmd)

## 注释
### 使用annovar
def annovar_anno(resultsDir, sampleID, threads, vcfFile):
	mkdir(resultsDir + "/annotation")
	humandb = "/path/to/humandb"
	buildver = "hg19"
	cmd = """
		convert2annovar.pl -format vcf4 {vcfFile} > {resultsDir}/tempFile/{sampleID}.avinput
		table_annovar.pl {resultsDir}/tempFile/{sampleID}.avinput \\
			{humandb} -buildver {buildver} \\
			-out {resultsDir}/annotation/{sampleID} -remove \\
			-protocol refGene,ensGene,cytoBand,avsnp150,gnomad211_genome,clinvar_20210207,cosmic_v92,dbnsfp41c \\
			-operation g,g,r,f,f,f,f,f \\
			-nastring . -thread {threads} -otherinfo
	""".format(vcfFile=vcfFile, resultsDir=resultsDir, sampleID=sampleID, humandb=humandb, threads=threads)
	os.system(cmd)

## CNV检测
### 使用cnvkit
### 使用cnvkit，需要先使用若干正常样本建立基线
#### 获得target文件
def cnvkit_make_target(bedFile, targetFile):
	refFlat = "refFlat.txt"
	cmd = """
		python cnvkit.py target {bedFile} \\
			--annotate {refFlat} \\
			-o {targetFile}
	""".format(bedFile=bedFile, refFlat=refFlat, targetFile=targetFile)
	os.system(cmd)

#### 获得antitarget文件
def cnvkit_make_antitarget(targetFile, accessTarget, antiTargetFile):
	reference = "ucsc.hg19.fasta"
	antiBedFile = "wgEncodeDukeMapabilityRegionsExcludable.bed"
	cmd = """
		python cnvkit.py access {reference} \\
			-x {antiBedFile} \\
			-o {accessTarget}
		python cnvkit.py antitarget \\
			{targetFile} -g {accessTarget} \\
			-o {antiTargetFile}
	""".format(reference=reference, antiBedFile=antiBedFile, accessTarget=accessTarget, antiTargetFile=antiTargetFile)
	os.system(cmd)

#### 使用正常样本建立基线
def cnvkit_make_reference(normalBam, targetFile, accessTarget, antiTargetFile, referenceCnn):
	refFlat = "refFlat.txt"
	reference = "ucsc.hg19.fasta"
	cmd = """
		python cnvkit.py batch \\
			-n {normalBam} \\
			--output-reference {referenceCnn} \\
			-t {targetFile} \\
			-a {antiTargetFile} \\
			-f {reference} \\
			-g {accessTarget}
	""".format(normalBam=normalBam, referenceCnn=referenceCnn, targetFile=targetFile, antiTargetFile=antiTargetFile, reference=reference, accessTarget=accessTarget)
	os.system(cmd)

#### 分析CNV
def cnvkit_call(resultsDir, sampleID, referenceCnn):
	mkdir(resultsDir + "/cnv")
	refFlat = "refFlat.txt"
	cmd = """
		python cnvkit.py batch \\
			{resultsDir}/bam/{sampleID}_bqsr.bam \\
			-r {referenceCnn} \\
			-d {resultsDir}/cnv/{sampleID} \\
			-m hybrid \\
			--annotate {refFlat}
	""".format(resultsDir=resultsDir, sampleID=sampleID, referenceCnn=referenceCnn, refFlat=refFlat)
	os.system(cmd)

## SV
### 使用lumpy进行检测
def lumpy_sv(resultsDir, sampleID, threads):
	mkdir(resultsDir + "/fusion")
	cmd = """
		samtools view -bh -F 1294 {resultsDir}/bam/{sampleID}.bam \\
			| samtools sort -@ {threads} - \\
			-o {resultsDir}/tempFile/{sampleID}.discordants.bam
		samtools view -h {resultsDir}/bam/{sampleID}.bam \\
			| extractSplitReads_BwaMem \\
			-i stdin \\
			| samtools view -bSh - \\
			| samtools sort -@ {threads} - \\
			-o {resultsDir}/tempFile/{sampleID}.splitters.bam
		lumpyexpress -B {resultsDir}/bam/{sampleID}.bam \\
			-D {resultsDir}/tempFile/{sampleID}.discordants.bam \\
			-S {resultsDir}/tempFile/{sampleID}.splitters.bam
			-O {resultsDir}/fusion/{sampleID}.lumpy.vcf
		svtyper-sso \\
			-i {resultsDir}/fusion/{sampleID}.lumpy.vcf \\
			-B {resultsDir}/bam/{sampleID}.bam \\
			--cores {threads} \\
			-o {resultsDir}/fusion/{sampleID}.gt.vcf
	""".format(resultsDir=resultsDir, sampleID=sampleID, threads=threads)
	os.system(cmd)

## MSI
### 使用msisensor2进行检测
def msisensor2_msi(resultsDir, sampleID):
	mkdir(resultsDir + "/msi")
	model = "models_hg19"
	cmd = """
		msisensor2 msi -M {model} -t {resultsDir}/bam/{sampleID}.bam -o {resultsDir}/msi
	""".format(model=model, resultsDir=resultsDir, sampleID=sampleID)
	os.system(cmd)

## TMB
### 使用自建脚本
def TMB_calculator(resultsDir, sampleID):
	annovarFile = open(resultsDir + "/annotation/" + sampleID + ".txt", "r")
	for line in annovarFile:
		if line.startswith("#"):
			header = line.split("\n")[0].split("\t")
		else:
			continue

	break
			

