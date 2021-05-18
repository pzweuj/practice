version 1.0
# pzw
# 20210517
# Mapping

# bwa
# https://github.com/lh3/bwa
# sambamba
# https://lomereiter.github.io/sambamba/
task Bwa {
    input {
        String sample
        File cleanRead1
        File cleanRead2
        Int threads
    }

    File reference = "/home/bioinfo/ubuntu/database/hg19/ucsc.hg19.fasta"
    File ref_dict = "/home/bioinfo/ubuntu/database/hg19/ucsc.hg19.dict"
    File ref_amb = "/home/bioinfo/ubuntu/database/hg19/ucsc.hg19.fasta.amb"
    File ref_ann = "/home/bioinfo/ubuntu/database/hg19/ucsc.hg19.fasta.ann"
    File ref_bwt = "/home/bioinfo/ubuntu/database/hg19/ucsc.hg19.fasta.bwt"
    File ref_fai = "/home/bioinfo/ubuntu/database/hg19/ucsc.hg19.fasta.fai"
    File ref_pac = "/home/bioinfo/ubuntu/database/hg19/ucsc.hg19.fasta.pac"
    File ref_sa = "/home/bioinfo/ubuntu/database/hg19/ucsc.hg19.fasta.sa"

    command <<<
        bwa mem -t ~{threads} \
            -R "@RG\tPL:illumina\tSM:~{sample}\tID:~{sample}" \
            -M ~{reference} ~{cleanRead1} ~{cleanRead2} \
            | sambamba view -f bam -t ~{threads} \
            -S /dev/stdin > ~{sample}.bam
        sambamba sort ~{sample}.bam \
            -t ~{threads} -o ~{sample}.sort.bam

    >>>

    output {
        File sortBam = "~{sample}.sort.bam"
        File sortBamBai = "~{sample}.sort.bam.bai"
    }

    runtime {
        docker: "pzweuj/mapping"
    }
}


# gencore
# https://github.com/OpenGene/gencore
# 当使用fastp识别UMI时使用此方法进行去重
task Gencore {
    input {
        String sample
        File sortBam
        File sortBamBai
    }

    File reference = "/home/bioinfo/ubuntu/database/hg19/ucsc.hg19.fasta"
    File ref_fai = "/home/bioinfo/ubuntu/database/hg19/ucsc.hg19.fasta.fai"
    File GENCORE = "/home/bioinfo/ubuntu/software/gencore/gencore"
    File SAMBAMBA = "/home/bioinfo/ubuntu/software/sambamba-0.8.0/sambamba"

    command <<<
        ~{GENCORE} -i ~{sortBam} \
            -r ~{reference} \
            -o ~{sample}.umi.bam \
            -u UMI -s 2 -d 1 \
            -j ~{sample}.json -h ~{sample}.html
        ~{SAMBAMBA} sort \
            -t {threads} \
            ~{sample}.umi.bam \
            -o ~{sample}.marked.bam -p
    >>>

    output {
        File jsonReport = "~{sample}.json"
        File htmlReport = "~{sample}.html"
        File umiBam = "~{sample}.umi.bam"
        File markBam = "~{sample}.marked.bam"
        File markBamBai = "~{sample}.marked.bam.bai"
    }
}


# sambamba
# https://lomereiter.github.io/sambamba/
task MarkDuplicates {
    input {
        String sample
        File sortBam
        File sortBamBai
        Int threads
    }

    File SAMBAMBA = "/home/bioinfo/ubuntu/software/sambamba-0.8.0/sambamba"

    command <<<
        ~{SAMBAMBA} markdup \
            ~{sortBam} \
            ~{sample}.marked.bam \
            -p --overflow-list-size 600000 \
            -t ~{threads}
    >>>

    output {
        File markBam = "~{sample}.marked.bam"
        File markBamBai = "~{sample}.marked.bam.bai"
    }
}


# GATK Recalibrator
# https://gatk.broadinstitute.org/hc/en-us
task Recalibrator {
    input {
        String sample
        File bam
        File bai
    }

    File reference = "/home/bioinfo/ubuntu/database/hg19/ucsc.hg19.fasta"
    File ref_fai = "/home/bioinfo/ubuntu/database/hg19/ucsc.hg19.fasta.fai"
    File ref_dict = "/home/bioinfo/ubuntu/database/hg19/ucsc.hg19.dict"
    File millsIndel = "/home/bioinfo/ubuntu/database/hg19/Mills_and_1000G_gold_standard.indels.hg19.sites.vcf"
    File millsIndelIdx = "/home/bioinfo/ubuntu/database/hg19/Mills_and_1000G_gold_standard.indels.hg19.sites.vcf.idx"
    File genome1k = "/home/bioinfo/ubuntu/database/hg19/1000G_phase1.indels.hg19.sites.vcf"
    File genome1kIdx = "/home/bioinfo/ubuntu/database/hg19/1000G_phase1.indels.hg19.sites.vcf.idx"
    File dbsnp = "/home/bioinfo/ubuntu/database/hg19/dbsnp_138.hg19.vcf"
    File dbsnpIdx = "/home/bioinfo/ubuntu/database/hg19/dbsnp_138.hg19.vcf.idx"

    command <<<
        gatk BaseRecalibrator \
            --known-sites ~{millsIndel} \
            --known-sites ~{genome1k} \
            --known-sites ~{dbsnp} \
            -R ~{reference} \
            -I ~{bam} \
            -O ~{sample}.recal.table
    >>>

    output {
        File table = "~{sample}.recal.table"
    }

    runtime {
        docker: "aperdriau/gatk4.2.0"
    }

}

# GATK Recalibrator
# https://gatk.broadinstitute.org/hc/en-us
task ApplyBQSR {
    input {
        String sample
        File bam
        File bai
        File table
    }

    File reference = "/home/bioinfo/ubuntu/database/hg19/ucsc.hg19.fasta"
    File ref_fai = "/home/bioinfo/ubuntu/database/hg19/ucsc.hg19.fasta.fai"
    File ref_dict = "/home/bioinfo/ubuntu/database/hg19/ucsc.hg19.dict"

    command <<<
        gatk ApplyBQSR \
            -R ~{reference} \
            --bqsr-recal-file ~{table} \
            -I ~{bam} \
            -O ~{sample}.Realign.bam
    >>>

    output {
        File realignBam = "~{sample}.Realign.bam"
        File realignBamBai = "~{sample}.Realign.bai"
    }

    runtime {
        docker: "aperdriau/gatk4.2.0"
    }
}


# bamdst
# https://github.com/shiquan/bamdst
# 用于统计bam比对效果
task Bamdst {
    input {
        String sample
        File bam
        File bai
        File bed
    }

    File BAMDST = "/home/bioinfo/ubuntu/software/bamdst/bamdst"

    command <<<
        mkdir ~{sample}_tmp
        ~{BAMDST} \
            -p ~{bed} \
            -o ~{sample}_tmp \
            ~{bam}
        
        python3 <<CODE
        bamdstReportFile = open("~{sample}_tmp/coverage.report", "r")
        bamdstReport = open("~{sample}.bamdst.txt", "w")
        bamdstReport.write("~{sample} Bamdst QC Report\n")
        for line in bamdstReportFile:
            if line.startswith("#"):
                continue
            else:
                lines = line.lstrip()
                bamdstReport.write(lines)
        bamdstReport.close()
        bamdstReportFile.close()
        CODE
    >>>

    output {
        File bamdstReport = "~{sample}.bamdst.txt"
    }
}