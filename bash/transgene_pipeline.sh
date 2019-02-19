# 20190219
# pzw
###########################
sample=H04728D
reference=CriGri.embl.fa
plasmid=CED1108
threads=8
cover=30
showbases=10
countFilter=20
gap=100
###########################
# :<<'BLOCK

echo "merge reference and plasmid"
cat reference/$reference reference/${plasmid}.fasta > reference/merge.fa

echo "make merge index"
time bwa index reference/merge.fa

# BLOCK'

mkdir cleandata
echo "filter"
time fastp -i rawdata/${sample}_1.fq.gz -I rawdata/${sample}_2.fq.gz \
	-o cleandata/${sample}.clean_1.fastq.gz -O cleandata/${sample}.clean_2.fastq.gz \
	-j cleandata/${sample}.json -h cleandata/${sample}.html -w $threads

mkdir bam
echo "mapping"
time bwa mem -t $threads -Y reference/merge.fa \
	cleandata/${sample}.clean_1.fastq.gz \
	cleandata/${sample}.clean_2.fastq.gz > bam/${sample}.sam

echo "transform to bam"
time samtools view -bSh bam/${sample}.sam > bam/${sample}.bam
rm bam/${sample}.sam

echo "sort bam file"
time samtools sort -@ $threads bam/${sample}.bam -o bam/${sample}.sorted.bam
rm bam/${sample}.bam
samtools index bam/${sample}.sorted.bam

echo "extract header"
time samtools view bam/${sample}.sorted.bam -h | grep "@" > bam/${sample}.header

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
time python transgene_analysis_v1.5.py bam/${sample}.${plasmid}.final.bam \
	results/${sample}.results.txt \
	$cover $showbases

echo "analysis junction pairs"
time python transgene_pair_v0.5.py results/${sample}.results.txt \
	results/${sample}.results.pair.txt $countFilter $gap

echo "task done"
