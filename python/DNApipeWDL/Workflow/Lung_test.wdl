version 1.0
# pzw
# 20210517
# Lung Cancer Tumor-only mode(FFPE)

import "QC.wdl" as qc
import "Mapping.wdl" as mapping
import "Mutation.wdl" as mutation
import "SV.wdl" as sv

workflow Lung {
    input {
        String sample
        File bam
        File bai
        Int threads
    }

    File bed = "/home/bioinfo/ubuntu/database/hg19/cnv/tmerge.bed"
    Int depth = "1000"
    Float maf = "0.02"

    call sv.Lumpy as Lumpy {
        input:
            sample = sample,
            threads = threads,
            bam = bam,
            bai = bai
    }

    call sv.LumpyFilter as LumpyFilter {
        input:
            sample = sample,
            bam = bam,
            bai = bai,
            vcf = Lumpy.vcf,
            depth = depth,
            maf = maf
    }


    call sv.Manta as Manta {
        input:
            sample = sample,
            threads = threads,
            bam = bam,
            bai = bai
    }


    call sv.MantaFilter as MantaFilter {
        input:
            sample = sample,
            vcf = Manta.vcf,
            depth = depth,
            maf = maf
    }

    call sv.FusionAnno as FusionAnnotation {
        input:
            sample = sample,
            pair = "",
            lumpyVcf = LumpyFilter.vcfFilter,
            mantaVcf = MantaFilter.vcfFilter
    }

}
