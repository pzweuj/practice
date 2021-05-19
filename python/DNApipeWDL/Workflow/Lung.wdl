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
        File rawRead1
        File rawRead2
        Int threads
    }

    File bed = "/home/bioinfo/ubuntu/database/hg19/cnv/tmerge.bed"
    Int depth = "1000"
    Float maf = "0.02"

    call qc.Fastp as QC {
        input:
            sample = sample,
            threads = threads,
            rawRead1 = rawRead1,
            rawRead2 = rawRead2
    }

    call qc.FastpFormat as FastqStat {
        input:
            sample = sample,
            jsonReport = QC.jsonReport
    }

    call mapping.Bwa as Mapping {
        input:
            sample = sample,
            threads = threads,
            cleanRead1 = QC.cleanRead1,
            cleanRead2 = QC.cleanRead2
    }

    call mapping.MarkDuplicates as MarkDup {
        input:
            sample = sample,
            threads = threads,
            sortBam = Mapping.sortBam,
            sortBamBai = Mapping.sortBamBai
    }

    call mapping.Recalibrator as Recal {
        input:
            sample = sample,
            bam = MarkDup.markBam,
            bai = MarkDup.markBamBai
    }

    call mapping.ApplyBQSR as ApplyBQSR {
        input:
            sample = sample,
            bam = MarkDup.markBam,
            bai = MarkDup.markBamBai,
            table = Recal.table
    }

    call mapping.Bamdst as BamStat {
        input:
            sample = sample,
            bam = ApplyBQSR.realignBam,
            bai = ApplyBQSR.realignBamBai,
            bed = bed
    }

    call mutation.Mutect2 as Mutect2 {
        input:
            sample = sample,
            bam = ApplyBQSR.realignBam,
            bai = ApplyBQSR.realignBamBai,
            threads = threads,
            bed = bed
    }

    call mutation.Filter as VcfFilter {
        input:
            sample = sample,
            vcf = Mutect2.vcf,
            depth = depth
    }

    call sv.Lumpy as Lumpy {
        input:
            sample = sample,
            threads = threads,
            bam = Mapping.sortBam,
            bai = Mapping.sortBamBai
    }

    call sv.Manta as Manta {
        input:
            sample = sample,
            threads = threads,
            bam = Mapping.sortBam,
            bai = Mapping.sortBamBai
    }

    call sv.LumpyFilter as LumpyFilter {
        input:
            sample = sample,
            bam = Mapping.sortBam,
            bai = Mapping.sortBamBai,
            vcf = Lumpy.vcf,
            depth = depth,
            maf = maf
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
