##################
RGSM=KPGP-00245
RGID=KPGP-00245
RGPU=KPGP
RGLB=KPGP
sample1=../Rawdata/KPGP-00245_L001_R1.fq.gz
sample2=../Rawdata/KPGP-00245_L001_R2.fq.gz
#################
export PATH=$PATH:/home/zw/software/WXS

# dataset
reference=/home/zw/database/b37/human_g1k_v37_decoy.fasta
indel1=/home/zw/database/b37/1000G_omni2.5.b37.vcf
indel2=/home/zw/database/b37/1000G_phase1.indels.b37.vcf


# mapping
bwa mem -t 8 $reference $sample1 $sample2 > ../Temp/${RGSM}.sam

# sorting
samtools view -Sb ../Temp/${RGSM}.sam > ../Temp/${RGSM}.bam
samtools sort -@ 8 -O bam -o ../Temp/${RGSM}.sorted.bam ../Temp/${RGSM}.bam

# mark duplicates
gatk MarkDuplicates \
	-I ../Temp/${RGSM}.sorted.bam \
	-O ../Temp/${RGSM}.marked_dup.bam \
	-M ../Bam/${RGSM}_marked_dup.txt \
	--REMOVE_DUPLICATES false

# add head info
gatk AddOrReplaceReadGroups \
	-I ../Temp/${RGSM}.marked_dup.bam \
	-O ../Temp/${RGSM}.addhead.bam \
	--RGID $RGID \
	--RGLB $RGLB \
	--RGPL illumina \
	--RGPU $RGPU \
	--RGSM $RGSM

# make index
gatk BuildBamIndex \
	-I ../Temp/${RGSM}.addhead.bam


# indel region realign
gatk BaseRecalibrator \
	-I ../Temp/${RGSM}.addhead.bam \
	-R $reference \
	--known-sites $indel1 \
	--known-sites $indel2 \
	-O ../Temp/${RGSM}_recal_data.table

# commit 
gatk ApplyBQSR \
	-R $reference \
	-I ../Temp/${RGSM}.addhead.bam \
	--bqsr-recal-file ../Temp/${RGSM}_recal_data.table \
	-O ../Bam/${RGSM}.final.bam

echo "**task done**"
