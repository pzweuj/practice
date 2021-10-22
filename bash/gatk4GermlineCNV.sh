#!/bin/bash
# pzw
# 20211022

func() {
    echo "Usage:"
    echo "gatk4GermlineCNV.sh -i <bamPath> -o <outputPath> -b <bed> -r <ref.fa> -p <contig.ploidy>"
    echo "Description"
    echo "bamPath, bam and bai files directory"
    echo "outputPath, output directory"
    echo "bedFile, bed file without header"
    echo "ref.fa, reference genome in fasta format"
    echo "contig.ploidy, ploidy file"
    exit -1
}

while getopts "i:o:b:r:p:" OPT; do
    case $OPT in
        i) bamDir="$OPTARG";;
        o) output="$OPTARG";;
        b) bed="$OPTARG";;
        r) reference="$OPTARG";;
        p) ploidy="$OPTARG";;
        h) func;;
        ?) func;;
    esac
done

# reference=$1
# bed=$2
# output=$3
# bamDir=$4
# ploidy=$5

# collect raw counts / --bin-length 0 no binning will be performed
if [ ! -d $output ]; then
    mkdir $output;
fi

gatk PreprocessIntervals \
    -R $reference \
    -L $bed \
    --bin-length 0 \
    --padding 0 \
    -imr OVERLAPPING_ONLY \
    -O $output/target.prep.interval_list

# annotate GC content
# gatk AnnotateIntervals \
#     -L $bed \
#     -R $reference \
#     -imr OVERLAPPING_ONLY \
#     -O $output/target.annotated.tsv

# loop
if [ ! -d $output/counts ]; then
    mkdir $output/counts;
fi
for i in $bamDir/*.bam;
do
    bamName=$(basename "$i" .bam);
    echo $bamName;

    # collect reads per bin
    gatk CollectReadCounts \
        -L $output/target.prep.interval_list \
        -R $reference \
        -imr OVERLAPPING_ONLY \
        -I $i \
        --format TSV \
        -O $output/counts/$bamName.tsv;    
done

# 获得query
if [ -f $output/tmp ]; then
    rm $output/tmp;
fi
ls $output/counts/*.tsv | while read line;
do
    echo -n "-I $line " >> $output/tmp;
done
sampleQuery=$(cat $output/tmp)
rm $output/tmp

# FilterIntervals Option 测试失败
# gatk FilterIntervals \
#     -L $output/target.prep.interval_list \
#     --annotated-intervals $output/target.annotated.tsv \
#     $sampleQuery \
#     -imr OVERLAPPING_ONLY \
#     -O $output/cohort.gc.filtered.interval_list

# call autosomal and allosomal contig ploidy
# need python package gcnvkernel/ python3.4 or later
gatk DetermineGermlineContigPloidy \
    -L $output/target.prep.interval_list \
    --interval-merging-rule OVERLAPPING_ONLY \
    $sampleQuery \
    --contig-ploidy-priors $ploidy \
    --output $output/ploidy \
    --output-prefix ploidy \
    --verbosity DEBUG

# call cnv COHORT mode
gatk GermlineCNVCaller \
    --run-mode COHORT \
    -L $output/target.prep.interval_list \
    $sampleQuery \
    --contig-ploidy-calls $output/ploidy/ploidy-calls \
    --interval-merging-rule OVERLAPPING_ONLY \
    --output $output/cohort \
    --output-prefix cohort \
    --verbosity DEBUG

# consolidate sample results
if [ ! -d $output/output ]; then
    mkdir $output/output;
fi
if [ -f $output/sample.list ]; then
    rm $output/sample.list;
fi
echo -e "#SM\tSAMPLEID" > $output/sample.list
for ((i=0; i<$((`ls $bamDir | wc -l` / 2)); i++));
do
    echo "sample_${i}";
    gatk PostprocessGermlineCNVCalls \
        --model-shard-path $output/cohort/cohort-model \
        --calls-shard-path $output/cohort/cohort-calls \
        --sample-index $i \
        --contig-ploidy-calls $output/ploidy/ploidy-calls \
        --output-genotyped-intervals $output/output/sample_${i}_cohort.vcf.gz \
        --output-genotyped-segments $output/output/sample_${i}_segment.cohort.vcf.gz \
        --output-denoised-copy-ratios $output/output/sample_${i}_ratio.txt;

    # rename / name from bam SM tag
    sampleName=`zcat $output/output/sample_${i}_cohort.vcf.gz | grep "#CHROM" | cut -f 10`;
    echo -e $sampleName"\tSample_"$i >> $output/sample.list;
    mv $output/output/sample_${i}_cohort.vcf.gz $output/output/${sampleName}_cohort.vcf.gz;
    mv $output/output/sample_${i}_cohort.vcf.gz.tbi $output/output/${sampleName}_cohort.vcf.gz.tbi;
    mv $output/output/sample_${i}_segment.cohort.vcf.gz $output/output/${sampleName}_segment.cohort.vcf.gz;
    mv $output/output/sample_${i}_segment.cohort.vcf.gz.tbi $output/output/${sampleName}_segment.cohort.vcf.gz.tbi;
    mv $output/output/sample_${i}_ratio.txt $output/output/${sampleName}_ratio.txt;
done
