#!/usr/bin/bash
######################################################
fq1=../Rawdata/NPC29F-N_1.fastq.gz
fq2=../Rawdata/NPC29F-N_2.fastq.gz
RGID=NPC29F
RGLB=NPC29F
RGSM=NPC29F-N
RGPU=unit1
######################################################
trimmomatic=/media/netdisk246/ZhaowenPan/WellPlay/software/Trimmomatic-0.36/trimmomatic-0.36.jar
gatk=/media/netdisk246/Bioinformatics/software/GenomeAnalysisTK/gatk-4.0.4.0/gatk

# 数据库
reference=/media/netdisk246/Bioinformatics/database/GATK/library/b37/human_g1k_v37_decoy.fasta
bundle=/media/netdisk246/Bioinformatics/database/GATK/library/b37

# 对数据库建立索引
$gatk IndexFeatureFile --feature-file $bundle/xxx.vcf

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
bwa mem -t 8 $reference $fq1 $fq2 > ../Temp/sample.bwa.sam

picard SortSam \
	I=../Temp/sample.bwa.sam \
	O=../Temp/sample.sorted.bam \
	SO=coordinate

# 标记重复序列
picard MarkDuplicates \
	I=../Temp/sample.sorted.bam \
	O=../Temp/sample.marked_dup.bam \
	M=../Temp/sample_marked_dup.txt \
	REMOVE_DUPLICATES=true

picard AddOrReplaceReadGroups \
	I=../Temp/sample.marked_dup.bam \
	O=../Temp/sample.addhead.bam \
	RGID=$RGID \
	RGLB=$RGLB \
	RGPL=illumina \
	RGPU=$RGPU \
	RGSM=$RGSM

picard BuildBamIndex \
	I=../Temp/sample.addhead.bam

# BQSR
time $gatk BaseRecalibrator \
	-R $reference \
	-I ../Temp/sample.addhead.bam \
	--known-sites $bundle/1000G_phase1.indels.b37.vcf \
	--known-sites $bundle/Mills_and_1000G_gold_standard.indels.b37.vcf \
	--known-sites $bundle/dbsnp_138.b37.vcf \
	-O ../Temp/sort.md.recal.table && \
	echo "** ${sample} BQSR done"

time $gatk ApplyBQSR \
	--bqsr-recal-file ../Temp/sort.md.recal.table \
	-R $reference \
	-I ../Temp/sample.addhead.bam \
	-O ../Bam/${RGSM}.final.bam && \
	echo "** Apply BQSR done **"

picard BuildBamIndex I=../Bam/${RGSM}.final.bam && \
	echo "** make pre bam done **"
