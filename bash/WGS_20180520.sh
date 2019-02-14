#!/usr/bin/bash
######################################################
fq1=../Rawdata/NPC29F-N_1.fastq.gz
fq2=../Rawdata/NPC29F-N_2.fastq.gz
RGID=NPC29F
RGLB=NPC29F
RGSM=NPC29F-N
RGPU=unit1
######################################################

#使用的软件
trimmomatic=/PATH/to/software/Trimmomatic-0.36/trimmomatic-0.36.jar
gatk=/PATH/to/software/GenomeAnalysisTK/gatk-4.0.4.0/gatk
bwa=/PATH/to/bwa
picard=/PATH/to/picard

# 数据库
reference=/PATH/to/GATK/library/b37/human_g1k_v37_decoy.fasta
bundle=/PATH/to/GATK/library/b37

# 对数据库建立索引
# $gatk IndexFeatureFile --feature-file $bundle/xxx.vcf

mkdir ../Temp
mkdir ../Bam

# # 使用Trimmomatic对原始数据进行质控
# time java -jar $trimmomatic PE \
# 	$fq1 $fq2 \
# 	../Temp/paired.1.fq.gz ../Temp/unpaired.1.fq.gz \
# 	../Temp/paired.2.fq.gz ../Temp/unpaired.2.fq.gz \
# 	ILLUMINACLIP:/media/netdisk246/ZhaowenPan/WellPlay/software/Trimmomatic-0.36/adapters/TruSeq3-PE-2.fa:2:30:10:8:True \
# 	SLIDINGWINDOW:5:15 LEADING:5 TRAILING:5 MINLEN:50 && echo "** fq QC done **"

# 使用bwa mem进行比对
# $bwa index $reference
$bwa mem -t 8 $reference $fq1 $fq2 > ../Temp/${RGSM}.bwa.sam

# 排序并且将sam文件转换为bam
java -jar $picard SortSam \
	I=../Temp/${RGSM}.bwa.sam \
	O=../Temp/${RGSM}.sorted.bam \
	SO=coordinate

# 标记重复序列
java -jar $picard MarkDuplicates \
	I=../Temp/${RGSM}.sorted.bam \
	O=../Temp/${RGSM}.marked_dup.bam \
	M=../Temp/${RGSM}_marked_dup.txt \
	REMOVE_DUPLICATES=true

# 添加头信息
java -jar $picard AddOrReplaceReadGroups \
	I=../Temp/${RGSM}.marked_dup.bam \
	O=../Temp/${RGSM}.addhead.bam \
	RGID=$RGID \
	RGLB=$RGLB \
	RGPL=illumina \
	RGPU=$RGPU \
	RGSM=$RGSM

# 创建索引
java -jar $picard BuildBamIndex \
	I=../Temp/${RGSM}.addhead.bam

# BQSR 质量值校正
time $gatk BaseRecalibrator \
	-R $reference \
	-I ../Temp/${RGSM}.addhead.bam \
	--known-sites $bundle/1000G_phase1.indels.b37.vcf \
	--known-sites $bundle/Mills_and_1000G_gold_standard.indels.b37.vcf \
	--known-sites $bundle/dbsnp_138.b37.vcf \
	-O ../Temp/${RGSM}.sort.md.recal.table && \
	echo "** ${RGSM} BQSR done"

# 运行校正
time $gatk ApplyBQSR \
	--bqsr-recal-file ../Temp/${RGSM}.sort.md.recal.table \
	-R $reference \
	-I ../Temp/${RGSM}.addhead.bam \
	-O ../Bam/${RGSM}.final.bam && \
	echo "** Apply BQSR done **"

# 创建索引
java -jar $picard BuildBamIndex I=../Bam/${RGSM}.final.bam && \
	echo "** make pre bam done **"

# Variants Calling
time $gatk HaplotypeCaller \
	-R $reference \
	-I ../Bam/${RGSM}.final.bam \
	-O ../Vcf/${RGSM}.HC.vcf.gz && \
	echo "** ${RGSM}.HC.vcf.gz done **"

# VQSR校正 SNP模式
time $gatk VariantRecalibrator \
	-R $reference \
	-V ../Vcf/${RGSM}.HC.vcf.gz \
	-resource hapmap,known=false,training=true,truth=true,prior=15.0:$bundle/hapmap_3.3_b37_pop_stratified_af.vcf.gz \
	-resource omni,known=false,training=true,truth=false,prior=12.0:$bundle/1000G_omni2.5.b37.vcf \
	-resource 1000G,known=false,training=true,truth=false,prior=10.0:$bundle/1000G_phase1.snps.high_confidence.b37.vcf \
	-resource dbsnp,known=true,training=false,truth=false,prior=6.0:$bundle/dbsnp_138.b37.vcf \
	-an DP -an QD -an FS -an SOR -an ReadPosRankSum -an MQRankSum \
	-mode SNP \
	-tranche 100.0 -tranche 99.9 -tranche 99.0 -tranche 95.0 -tranche 90.0 \
	--rscript-file ../Vcf/${RGSM}.HC.snps.plots.R \
	--tranches-file ../Temp/${RGSM}.HC.snps.tranches \
	-O ../Temp/${RGSM}.snp.recal

# 确认VQSR校正 SNP模式
time $gatk ApplyVQSR \
	-R $reference \
	-V ../Vcf/${RGSM}.HC.vcf.gz \
	-ts-filter-level 99.0 \
	--tranches-file ../Temp/${RGSM}.HC.snps.tranches \
	--recal-file ../Temp/${RGSM}.snp.recal \
	-mode SNP \
	-O ../Temp/${RGSM}.snp.vcf.gz

# VQSR校正 INDEL模式
time $gatk VariantRecalibrator \
	-R $reference \
	-V ../Temp/${RGSM}.snp.vcf.gz \
	-resource mills,known=true,training=true,truth=true,prior=12.0:$bundle/Mills_and_1000G_gold_standard.indels.b37.vcf \
	-an DP -an QD -an FS -an SOR -an ReadPosRankSum -an MQRankSum \
	-mode INDEL \
	--max-gaussians 6 \
	-rscript-file ../Vcf/${RGSM}.HC.indels.plots.R \
	--tranches-file ../Temp/${RGSM}.HC.indels.tranches \
	-O ../Temp/${RGSM}.indels.recal

# 确认VQSR校正 INDEL模式
time $gatk ApplyVQSR \
	-R $reference \
	-V ../Temp/${RGSM}.snp.vcf.gz \
	-ts-filter-level 99.0 \
	--tranches-file ../Temp/${RGSM}.HC.indels.tranches \
	--recal-file ../Temp/${RGSM}.indels.recal \
	-mode INDEL \
	-O ../Vcf/${RGSM}.VQSR.vcf.gz && \
	echo "VQSR done"

# 删除中间文件
# rm -rf ../Temp