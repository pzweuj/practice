# pzw
# 20180324

normal1=../Rawdata/NPC29F-N_1.fastq.gz
normal2=../Rawdata/NPC29F-N_2.fastq.gz
tumor1=../Rawdata/NPC29F-T_1.fastq.gz
tumor2=../Rawdata/NPC29F-T_2.fastq.gz

reference=/media/netdisk246/Bioinformatics/database/GATK/library/b37/human_g1k_v37_decoy.fasta
indel1=/media/netdisk246/Bioinformatics/database/GATK/library/b37/1000G_phase1.indels.b37.vcf
indel2=/media/netdisk246/Bioinformatics/database/GATK/library/b37/Mills_and_1000G_gold_standard.indels.b37.vcf
snp=/media/netdisk246/Bioinformatics/database/GATK/library/b37/1000G_phase1.snps.high_confidence.b37.vcf
af=/media/netdisk246/Bioinformatics/database/GATK/library/b37/af-only-gnomad.raw.sites.b37.vcf.gz

gatk=/media/netdisk246/Bioinformatics/software/GenomeAnalysisTK/gatk-4.0.2.1/gatk

# tumor bam
mkdir ../Temp

bwa mem -t 6 $reference $tumor1 $tumor2 > ../Temp/tumor.sam

picard SortSam \
	I=../Temp/tumor.sam \
	O=../Temp/tumor.sorted.bam \
	SO=coordinate

picard MarkDuplicates \
	I=../Temp/tumor.sorted.bam \
	O=../Temp/tumor.marked_dup.bam \
	M=../Temp/tumor_marked_dup.txt \
	REMOVE_DUPLICATES=true

picard AddOrReplaceReadGroups \
	I=../Temp/tumor.marked_dup.bam \
	O=../Temp/tumor.addhead.bam \
	RGID=NPC29F \
    RGLB=lib1 \
    RGPL=illumina \
    RGPU=unit1 \
    RGSM=tumor

picard BuildBamIndex \
	I=../Temp/tumor.addhead.bam

$gatk BaseRecalibrator \
	-I ../Temp/tumor.addhead.bam \
	-R $reference \
	--known-sites $indel1 \
	--known-sites $indel2 \
	-O ../Temp/tumor_recal_data.table

$gatk ApplyBQSR \
	-R $reference \
	-I ../Temp/tumor.addhead.bam \
	--bqsr-recal-file ../Temp/tumor_recal_data.table \
	-O ../Bam/tumor.final.bam

rm -rf ../Temp

# make normal bam
mkdir ../Temp

bwa mem -t 6 $reference $normal1 $normal2 > ../Temp/normal.sam

picard SortSam \
	I=../Temp/normal.sam \
	O=../Temp/normal.sorted.bam \
	SO=coordinate

picard MarkDuplicates \
	I=../Temp/normal.sorted.bam \
	O=../Temp/normal.marked_dup.bam \
	M=../Temp/normal_marked_dup.txt \
	REMOVE_DUPLICATES=true

picard AddOrReplaceReadGroups \
	I=../Temp/normal.marked_dup.bam \
	O=../Temp/normal.addhead.bam \
	RGID=NPC29F \
	RGLB=lib1 \
	RGPL=illumina \
	RGPU=unit1 \
	RGSM=normal

picard BuildBamIndex \
	I=../Temp/normal.addhead.bam

$gatk BaseRecalibrator \
	-I ../Temp/normal.addhead.bam \
	-R $reference \
	--known-sites $indel1 \
	--known-sites $indel2 \
	-O ../Temp/normal_recal_data.table

$gatk ApplyBQSR \
	-R $reference \
	-I ../Temp/normal.addhead.bam \
	--bqsr-recal-file ../Temp/normal_recal_data.table \
	-O ../Bam/normal.final.bam

rm -rf ../Temp

# mutect2
mkdir ../Temp

$gatk Mutect2 \
	-R $reference \
	-I ../Bam/normal.final.bam \
	-tumor normal \
	--germline-resource $af \
	-O ../Temp/normal_for_pon.vcf.gz

$gatk CreateSomaticPanelOfNormals \
	-vcfs ../Temp/normal_for_pon.vcf.gz \
	-O ../Temp/pon.vcf.gz

$gatk Mutect2 \
	-R $reference \
	-I ../Bam/tumor.final.bam \
	-I ../Bam/normal.final.bam \
	-tumor tumor \
	-normal normal \
	-pon ../Temp/pon.vcf.gz\
	--germline-resource $af \
	--af-of-alleles-not-in-resource 0.0000025 \
	--disable-read-filter MateOnSameContigOrNoMappedMateReadFilter \
	-O ../Vcf/somatic_m2.vcf.gz \
	-bamout ../Bam/tumor_normal_m2.bam

$gatk GetPileupSummaries \
	-I ../Bam/tumor.final.bam \
	-V $snp \
	-V $indel1 \
	-V $indel2 \
	-O ../Temp/summaries.table

$gatk CalculateContamination \
	-I ../Temp/summaries.table \
	-O ../Temp/calculatecontamination.table

$gatk FilterMutectCalls \
	-V ../Vcf/somatic_m2.vcf.gz \
	--contamination-table ../Temp/calculatecontamination.table \
	-O ../Vcf/somatic_oncefiltered.vcf.gz

$gatk CollectSequencingArtifactMetrics \
	-I ../Bam/tumor.final.bam \
	-O ../Temp/tumor_artifact \
	-EXT ".txt" \
	-R $reference

$gatk FilterByOrientationBias \
	-V ../Vcf/somatic_oncefiltered.vcf.gz \
	--artifact-modes 'G/T' \
	-P ../Temp/tumor_artifact.pre_adapter_detail_metrics.txt \
	-O ../Vcf/somatic_twicefiltered.vcf.gz

rm -rf ../Temp
