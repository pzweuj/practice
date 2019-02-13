sample=ERR900939

mkdir cleandata
fastp -i rawdata/${sample}_1.fastq.gz \
	-I rawdata/${sample}_2.fastq.gz \
	-o cleandata/${sample}.clean_1.fastq.gz \
	-O cleandata/${sample}.clean_2.fastq.gz

echo "filter done"

mkdir map2humanInflu
bwa mem -t 8 influenza.human.fa cleandata/${sample}.clean_1.fastq.gz cleandata/${sample}.clean_2.fastq.gz > map2humanInflu/${sample}.sam

echo "mapping done"

samtools view -bSh -F 12 map2humanInflu/${sample}.sam > map2humanInflu/${sample}.map.bam

echo "filter unmapped reads"

bedtools bamtofastq -i map2humanInflu/${sample}.map.bam -fq map2humanInflu/${sample}.map.fastq

echo "trans to fastq"

cat map2humanInflu/${sample}.map.fastq | paste - - - - | sed 's/^@/>/g'| cut -f1-2 | tr '\t' '\n' > map2humanInflu/${sample}.map.fasta

echo "trans to fasta"

mkdir blastResults
blastn -query map2humanInflu//${sample}.map.fasta \
        -db influenza_select/influenza.human.HA.select.fa \
        -out blastResults/${sample}-results-HA.txt \
        -evalue 10

blastn -query map2humanInflu//${sample}.map.fasta \
        -db influenza_select/influenza.human.NA.select.fa \
        -out blastResults/${sample}-results-NA.txt \
        -evalue 10


echo "blast done"

python get_influenza_subtype.py blastResults/${sample}-results-HA.txt blastResults/${sample}-results-NA.txt blastResults/${sample}-results.txt

echo "please see results"






