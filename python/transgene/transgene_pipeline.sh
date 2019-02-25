# 20190218
# pzw
###########################
sample=H04726D
reference=CriGri.embl.fa
plasmid=CED1108
threads=8
cover=30
showbases=10
countFilter=0
gap=100
###########################
# :<<'BLOCK

# echo "merge reference and plasmid"
# cat reference/$reference reference/${plasmid}.fa > reference/merge.fa

# echo "make merge index"
# time bwa index reference/merge.fa

# BLOCK'

mkdir rawdata
ln -s /oridata/origin-data-backup/Internal-OriData/006-NuoHe-Novaseq-QB058-20190221/Rawdata/QB058-${sample}T* \
	rawdata/

zcat rawdata/QB058-${sample}T5_L1_1.fq.gz rawdata/QB058-${sample}T5_L2_1.fq.gz rawdata/QB058-${sample}T5_L3_1.fq.gz \
	rawdata/QB058-${sample}T7_L1_1.fq.gz rawdata/QB058-${sample}T7_L2_1.fq.gz rawdata/QB058-${sample}T7_L3_1.fq.gz | \
	gzip - > rawdata/${sample}_1.fq.gz

zcat rawdata/QB058-${sample}T5_L1_2.fq.gz rawdata/QB058-${sample}T5_L2_2.fq.gz rawdata/QB058-${sample}T5_L3_2.fq.gz \
	rawdata/QB058-${sample}T7_L1_2.fq.gz rawdata/QB058-${sample}T7_L2_2.fq.gz rawdata/QB058-${sample}T7_L3_2.fq.gz | \
	gzip - > rawdata/${sample}_2.fq.gz


mkdir cleandata
echo "filter"
time fastp -i rawdata/${sample}_1.fq.gz -I rawdata/${sample}_2.fq.gz \
	-o cleandata/${sample}.clean_1.fq.gz -O cleandata/${sample}.clean_2.fq.gz \
	-j cleandata/${sample}.json -h cleandata/${sample}.html -w $threads

mkdir bam
echo "mapping"
time bwa mem -t $threads -Y ~/workspace/database/CriGri/CriGri_${plasmid}.fasta \
	cleandata/${sample}.clean_1.fq.gz \
	cleandata/${sample}.clean_2.fq.gz > bam/${sample}.sam

echo "transform to bam"
time samtools view -bSh bam/${sample}.sam > bam/${sample}.bam
rm bam/${sample}.sam

echo "sort bam file"
time samtools sort -@ $threads bam/${sample}.bam -o bam/${sample}.sorted.bam
rm bam/${sample}.bam
samtools index bam/${sample}.sorted.bam

echo "extract header"
time samtools view bam/${sample}.sorted.bam -h | grep "@SQ" > bam/${sample}.header

echo "extract junction reads"
time samtools view bam/${sample}.sorted.bam | grep "SA" | grep $plasmid > bam/${sample}.${plasmid}.sam

echo "merge reads and header"
cat bam/${sample}.header bam/${sample}.${plasmid}.sam > bam/${sample}.${plasmid}.final.sam
rm bam/${sample}.header bam/${sample}.${plasmid}.sam

echo "transform to bam"
time samtools view -bSh bam/${sample}.${plasmid}.final.sam > bam/${sample}.${plasmid}.final.bam
rm bam/${sample}.${plasmid}.final.sam
samtools index bam/${sample}.${plasmid}.final.bam

mkdir results
echo "analysis junction reads"
time python ~/workspace/bin/transgene_analysis_v1.6.py bam/${sample}.${plasmid}.final.bam \
	results/${sample}.results.txt \
	$cover $showbases

echo "analysis junction pairs"
time python ~/workspace/bin/transgene_pair_v0.5.py results/${sample}.results.txt \
	results/${sample}.results.pair.txt $countFilter $gap

echo "annotate gene region"
time python ~/workspace/bin/transgene_annotate_v1.1.py results/${sample}.results.pair.txt \
	results/${sample}.results.anno.txt \
	/home/zhaowen/workspace/database/CriGri/CriGri_annoDB.txt

echo "task done"
