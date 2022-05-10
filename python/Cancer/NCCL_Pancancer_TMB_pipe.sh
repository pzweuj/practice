#!/usr/bin/bash
sample=$1

fastp -i Fastq/${sample}.R1.fastq.gz \
    -I Fastq/${sample}.R2.fastq.gz \
    -o cleandata/${sample}_R1.fq.gz \
    -O cleandata/${sample}_R2.fq.gz \
    -j QC/${sample}.json \
    -h QC/${sample}.html \
    -w 16
bwa mem -t 16 -R "@RG\tPL:illumina\tSM:${sample}\tPU:YK\tID:${sample}" \
    -v 1 -M ucsc.hg19.fasta \
    cleandata/${sample}_R1.fq.gz cleandata/${sample}_R2.fq.gz \
    > bam/${sample}.sam
sambamba view bam/${sample}.sam -S -h \
    -f bam -p -t 16 -o bam/${sample}.bam
mkdir bam/${sample}_tmp
sambamba sort -t 16 bam/${sample}.bam -o bam/${sample}.sort.bam -p \
    --tmpdir bam/${sample}_tmp -m 64G
rm bam/${sample}.bam bam/${sample}.sam
gatk MarkDuplicates \
    -I bam/${sample}.sort.bam \
    -O bam/${sample}.marked.bam \
    -M bam/${sample}.marked.txt \
    --TMP_DIR bam/${sample}_tmp \
    --CREATE_INDEX true
mv bam/${sample}.marked.bai bam/${sample}.marked.bam.bai
gatk BaseRecalibrator \
    --known-sites Mills_and_1000G_gold_standard.indels.hg19.sites.vcf \
    --known-sites 1000G_phase1.indels.hg19.sites.vcf \
    --known-sites dbsnp_138.hg19.vcf \
    -R ucsc.hg19.fasta \
    -I bam/${sample}.marked.bam \
    -O bam/${sample}.recal.table
gatk ApplyBQSR \
    -R ucsc.hg19.fasta \
    --bqsr-recal-file bam/${sample}.recal.table \
    -I bam/${sample}.marked.bam \
    -O bam/${sample}.Realign.bam
mv bam/${sample}.Realign.bai bam/${sample}.Realign.bam.bai
rm bam/${sample}.marked.bam*
rm bam/${sample}.sort.bam*
rm bam/${sample}.marked.txt
rm -rf bam/${sample}_tmp

gatk GetPileupSummaries \
    -I bam/${sample}.Realign.bam \
    -V small_exac_common_3_hg19.vcf \
    -L small_exac_common_3_hg19.vcf \
    -O vcf/${sample}_tmp/${sample}.pipeups.table

mkdir vcf/${sample}_tmp/tmp

gatk CalculateContamination \
    -I vcf/${sample}_tmp/${sample}.pipeups.table \
    -matched vcf/2022NC.pipeups.table \
    -O vcf/${sample}_tmp/${sample}.contamination.table \
    --tumor-segmentation vcf/${sample}_tmp/${sample}.tumor_segments.txt \
    --tmp-dir vcf/${sample}_tmp/tmp

gatk Mutect2 \
    -R ucsc.hg19.fasta \
    -I bam/${sample}.Realign.bam \
    -I bam/2022NC.Realign.bam \
    -O vcf/${sample}.mutect2.vcf \
    -tumor ${sample} \
    -normal 2022NC \
    --germline-resource af-only-gnomad.raw.sites.hg19.vcf.gz \
    --native-pair-hmm-threads 8 \
    -L x.bed \
    -A Coverage -A GenotypeSummaries \
    -mbq 15 --force-active true \
    --callable-depth 150 \
    --f1r2-tar-gz vcf/${sample}_tmp/${sample}.f1r2.tar.gz \
    --max-reads-per-alignment-start 0

gatk LearnReadOrientationModel \
    -I vcf/${sample}_tmp/${sample}.f1r2.tar.gz \
    -O vcf/${sample}_tmp/${sample}.read-orientation-model.tar.gz \
    --tmp-dir vcf/${sample}_tmp/tmp

gatk FilterMutectCalls \
    -R ucsc.hg19.fasta \
    -V vcf/${sample}.mutect2.vcf \
    -L x.bed \
    --contamination-table vcf/${sample}_tmp/${sample}.contamination.table \
    --tumor-segmentation vcf/${sample}_tmp/${sample}.tumor_segments.txt \
    --orientation-bias-artifact-priors vcf/${sample}_tmp/${sample}.read-orientation-model.tar.gz \
    -O vcf/${sample}.mutect2.filter.vcf \
    --tmp-dir vcf/${sample}_tmp/tmp

bcftools view -s $sample vcf/${sample}.mutect2.filter.vcf > vcfFilter/${sample}.raw.vcf
gatk LeftAlignAndTrimVariants \
    -R ucsc.hg19.fasta \
    --split-multi-allelics \
    --keep-original-ac \
    -no-trim \
    -V vcfFilter/${sample}.raw.vcf -O vcfFilter/${sample}.left.vcf

java -jar ~/software/snpEff-5.0d/snpEff.jar \
    -c ~/software/snpEff-5.0d/snpEff.config -noStats hg19 \
    vcfFilter/${sample}.left.vcf > annotations/${sample}.snpeff.vcf

convert2annovar.pl -format vcf4 -allsample -withfreq \
    annotations/${sample}.snpeff.vcf \
    --includeinfo > annotations/${sample}.avinput

table_annovar.pl annotations/${sample}.avinput \
    humandb -buildver hg19 \
    -out annotations/${sample} -remove \
    -protocol refGene,avsnp150,gnomad211_genome,1000g2015aug_all,exac03,clinvar_20220320,JaxCkb,Civic,OncoKB,dbnsfp42a,dbscsnv11,cosmic92_coding,intervar_20180118,SpliceAI \
    -operation g,f,f,f,f,f,f,f,f,f,f,f,f,f \
    -nastring - -thread 8 -otherinfo
rm annotations/${sample}.avinput

python3 annoFixNCCL.py annotations/${sample}.hg19_multianno.txt annotations/${sample}.anno.txt refTranscript.txt

# end

