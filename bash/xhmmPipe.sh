#!/bin/bash
# 20211020
# 原始输入是经比对、排序、去重、BQSR的bam

echo "xhmmPipe.sh -h"
func() {
    echo "Usage:"
    echo "xhmmPipe.sh -i <bamPath> -o <outputPath> -b <bedFile> -r <ref.fa> -s <seqdb>"
    echo "Description"
    echo "bamPath, bam and bai files directory"
    echo "outputPath, output directory"
    echo "bedFile, bed file without header"
    echo "ref.fa, reference genome in fasta format"
    echo "seqdb, seqdb file"
    exit -1
}
while getopts "i:o:b:r:s:" OPT; do
    case $OPT in
        i) bamPath="$OPTARG";;
        o) output="$OPTARG";;
        b) bed="$OPTARG";;
        r) reference="$OPTARG";;
        s) seqdb="$OPTARG";;
        h) func;;
        ?) func;;
    esac
done

###############################
# bamPath=$1
# bed=$2
# reference=$3
# output=$4
# seqdb=$5
tmpdir=$output/tmp
###############################


# 创建tmpdir
if [ ! -d ${output} ]; then
    mkdir ${output}
fi
if [ ! -d ${tmpdir} ]; then
    mkdir ${tmpdir}
fi
if [ -f ${tmpdir}/bams.list ]; then
    rm ${tmpdir}/bams.list
fi
ls $bamPath | grep -v '.bai' | while read line; do echo ${bamPath}/${line} >> ${tmpdir}/bams.list; done

# 使用GATK3统计深度
java -jar /software/gatk3.8-1/GenomeAnalysisTK.jar \
    -T DepthOfCoverage -I ${tmpdir}/bams.list -L ${bed} \
    -R ${reference} \
    -dt BY_SAMPLE -dcov 5000 -l INFO --omitDepthOutputAtEachBase --omitLocusTable \
    --minBaseQuality 0 --minMappingQuality 20 --start 1 --stop 5000 --nBins 200 \
    --includeRefNSites \
    --countType COUNT_FRAGMENTS \
    -o ${tmpdir}/group

# 合并
xhmm --mergeGATKdepths \
    -o ${output}/data.rd.txt \
    --GATKdepths ${tmpdir}/group.sample_interval_summary

# 获得target区域GC比例
java -jar /software/gatk3.8-1/GenomeAnalysisTK.jar \
    -T GCContentByInterval -L ${bed} \
    -R ${reference} \
    -o ${tmpdir}/data.locus.gc.txt
cat ${tmpdir}/data.locus.gc.txt | awk '{if ($2 < 0.1 || $2 > 0.9) print $1}' > ${tmpdir}/extreme.gc.target.txt

# 可选，计算重复碱基比例并过滤，部分shell环境会运行失败
# java -jar /software/gatk3.8-1/picard.jar BedToIntervalList \
#     I=${bed} \
#     O=${tmpdir}/tmp.interval_list \
#     SD=$(echo ${reference} | sed 's/.fa/ /' | cut -d ' ' -f 1)".dict"
cat ${bed} | awk '{print $1":"$2"-"$3}' > ${tmpdir}/tmp.interval_list

interval_list_to_pseq_reg \
    ${tmpdir}/tmp.interval_list > ${tmpdir}/tmp.target.reg

pseq . loc-load \
    --locdb ${tmpdir}/tmp.target.locdb \
    --file ${tmpdir}/tmp.target.reg \
    --group targets \
    --out ${tmpdir}/tmp.target.locdb.loc-load \
    --noweb

pseq . loc-stats \
    --locdb ${tmpdir}/tmp.target.locdb \
    --group targets --seqdb ${seqdb} | \
    awk '{if (NR > 1) print $_}' | \
    sort -k1 -g | awk '{print $10}' | \
    paste ${tmpdir}/tmp.interval_list - | \
    awk '{print $1"\t"$2}' \
    > ${tmpdir}/data.locus_complexity.txt

cat ${tmpdir}/data.locus_complexity.txt | awk '{if ($2 > 0.25) print $1}' > ${tmpdir}/low_complexity_targets.txt

# 过滤求均值
xhmm --matrix \
    -r ${output}/data.rd.txt --centerData --centerType target \
    -o ${tmpdir}/data.filtered_centered.RD.txt \
    --outputExcludedTargets ${tmpdir}/data.filtered_centered.RD.txt.filtered_targets.txt \
    --outputExcludedSamples ${tmpdir}/data.filtered_centered.RD.txt.filtered_samples.txt \
    --excludeTargets ${tmpdir}/extreme.gc.target.txt --excludeTargets ${tmpdir}/low_complexity_targets.txt \
    --minTargetSize 10 --maxTargetSize 10000 \
    --minMeanTargetRD 10 --maxMeanTargetRD 500 \
    --minMeanSampleRD 25 --maxMeanSampleRD 200 \
    --maxSdSampleRD 150

# PCA聚类分析
xhmm --PCA -r ${tmpdir}/data.filtered_centered.RD.txt --PCAfiles ${output}/data.RD.PCA

# 降噪
xhmm --normalize -r ${tmpdir}/data.filtered_centered.RD.txt --PCAfiles ${output}/data.RD.PCA \
    --normalizeOutput ${output}/data.PCA_normalized.txt \
    --PCnormalizeMethod PVE_mean --PVE_mean_factor 0.7

# 过滤 Z检验
xhmm --matrix \
    -r ${output}/data.PCA_normalized.txt --centerData --centerType sample --zScoreData \
    -o ${output}/data.PCA_normalized.filtered.sample_zscores.RD.txt \
    --outputExcludedTargets ${tmpdir}/data.PCA_normalized.filtered.sample_zscores.RD.txt.filtered_targets.txt \
    --outputExcludedSamples ${tmpdir}/data.PCA_normalized.filtered.sample_zscores.RD.txt.filtered_samples.txt \
    --maxSdTargetRD 30

# 过滤原始深度
xhmm --matrix \
    -r ${output}/data.rd.txt \
    --excludeTargets ${tmpdir}/data.filtered_centered.RD.txt.filtered_targets.txt \
    --excludeTargets ${tmpdir}/data.PCA_normalized.filtered.sample_zscores.RD.txt.filtered_targets.txt \
    --excludeSamples ${tmpdir}/data.filtered_centered.RD.txt.filtered_samples.txt \
    --excludeSamples ${tmpdir}/data.PCA_normalized.filtered.sample_zscores.RD.txt.filtered_samples.txt \
    -o ${output}/data.same_filtered.RD.txt

# CNV检测
xhmm --discover \
    -p /software/statgen-xhmm-998f7c405974/params.txt \
    -r ${output}/data.PCA_normalized.filtered.sample_zscores.RD.txt -R ${output}/data.same_filtered.RD.txt \
    -c ${output}/data.xcnv -a ${output}/data.aux_xcnv -s ${output}/data

# 注释
xhmm --genotype \
    -p /software/statgen-xhmm-998f7c405974/params.txt \
    -r ${output}/data.PCA_normalized.filtered.sample_zscores.RD.txt -R ${output}/data.same_filtered.RD.txt \
    -g ${output}/data.xcnv -F ${reference} \
    -v ${output}/data.xhmm.vcf
