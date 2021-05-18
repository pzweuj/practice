version 1.0
# pzw
# 20210518
# SNV InDel


# Mutect2
# https://github.com/broadinstitute/gatk
task Mutect2 {
    input {
        File bam
        File bai
        String sample
        Int threads
        File? bed
    }

    File reference = "/home/bioinfo/ubuntu/database/hg19/ucsc.hg19.fasta"
    File ref_dict = "/home/bioinfo/ubuntu/database/hg19/ucsc.hg19.dict"
    File ref_fai = "/home/bioinfo/ubuntu/database/hg19/ucsc.hg19.fasta.fai"
    File gnomad = "/home/bioinfo/ubuntu/database/hg19/af-only-gnomad.raw.sites.hg19.vcf.gz"
    File gnomadIdx = "/home/bioinfo/ubuntu/database/hg19/af-only-gnomad.raw.sites.hg19.vcf.gz.tbi"

    command <<<
        gatk Mutect2 \
            -R ~{reference} \
            -I ~{bam} \
            -O ~{sample}.mutect2.vcf \
            -tumor ~{sample} \
            --germline-resource ~{gnomad} \
            --native-pair-hmm-threads ~{threads} \
            -L ~{bed} \
            -A Coverage -A GenotypeSummaries \
            --genotype-germline-sites false \
            --max-reads-per-alignment-start 0
    >>>

    output {
        File vcf = "~{sample}.mutect2.vcf"
    }

    runtime {
        docker: "aperdriau/gatk4.2.0"
    }
}


task HaplotypeCaller {
    input {
        String sample
        File bam
        File bai
        Int threads
        File? bed
    }

    File reference = "/home/bioinfo/ubuntu/database/hg19/ucsc.hg19.fasta"
    File ref_dict = "/home/bioinfo/ubuntu/database/hg19/ucsc.hg19.dict"

    command <<<
        gatk HaplotypeCaller \
            -R ~{reference} \
            -I ~{bam} \
            -O ~{sample}.hap.vcf \
            --native-pair-hmm-threads ~{threads} \
            -L ~{bed}
    >>>

    output {
        File vcf = "~{sample}.hap.vcf"
    }

    runtime {
        docker: "aperdriau/gatk4.2.0"
    }

}


# VarScan2
# http://varscan.sourceforge.net/
task VarScan2 {
    input {
        String sample
        String pair
        File tumorBam
        File tumorBai
        File normalBam
        File normalBai
        File bed
        Int depth
        Float maf
    }

    File VARSCAN = "/home/bioinfo/ubuntu/software/VarScan.v2.3.9/VarScan.v2.3.9.jar"
    File SAMTOOLS = "/home/bioinfo/ubuntu/software/samtools-1.11/samtools"
    File BCFTOOLS = "/home/bioinfo/ubuntu/software/bcftools-1.11/bcftools"
    File reference = "/home/bioinfo/ubuntu/database/hg19/ucsc.hg19.fasta"
    File ref_dict = "/home/bioinfo/ubuntu/database/hg19/ucsc.hg19.dict"
    File ref_fai = "/home/bioinfo/ubuntu/database/hg19/ucsc.hg19.fasta.fai"
    Int minCov = depth * maf

    command <<<
        ~{SAMTOOLS} mpileup -B -f ~{reference} -q 15 -d 10000 \
            ~{normalBam} ~{tumorBam} \
            | java -jar ~{VARSCAN} somatic -mpileup \
            ~{sample} \
            --min-coverage-normal ~{minCov} \
            --min-coverage-tumor ~{minCov} \
            --min-var-freq ~{maf} --strand-filter 1 --output-vcf
        ~{BCFTOOLS} reheader -f ~{ref_fai} ~{sample}.indel.vcf -o ~{sample}.indel.varscan.vcf
        ~{BCFTOOLS} reheader -f ~{ref_fai} ~{sample}.snp.vcf -o ~{sample}.snp.varscan.vcf
    >>>

    output {
        File indelVcf = "~{sample}.indel.varscan.vcf"
        File snpVcf = "~{sample}.snp.varscan.vcf"
    }
}

task VarScan2Filter {
    input {
        String sample
        String pair
        File snpVcf
        File indelVcf
    }

    File BCFTOOLS = "/home/bioinfo/ubuntu/software/bcftools-1.11/bcftools"

    command <<<
        python3 <<CODE
        import os
        outputFile = "~{sample}.merge.vcf"
        output = open(outputFile, "w")
        indel = open("~{indelVcf}", "r")
        snp = open("~{snpVcf}", "r")

        for i in indel:
            if i.startswith("#"):
                if "NORMAL\tTUMOR" in i:
                    i = i.replace("NORMAL\tTUMOR", "~{sample}\t~{pair}")
            else:
                ii = i.replace("\n", "").split("\t")
                li = ii[0:9]
                li.append(ii[10])
                li.append(ii[9])
                i = "\t".join(li) + "\n"
            output.write(i)
        indel.close()

        for s in snp:
            if not s.startswith("#"):
                ss = s.replace("\n", "").split("\t")
                si = ss[0:9]
                si.append(ss[10])
                si.append(ss[9])
                s = "\t".join(si) + "\n"
                output.write(s)
        snp.close()
        output.close()

        cmd = """
            {bcftools} sort {outputFile} -O v -o {sampleID}.vcf
        """.format(bcftools="~{BCFTOOLS}", outputFile=outputFile, sampleID="~{sample}", pairID="~{pair}")
        print(cmd)
        os.system(cmd)
        CODE
    >>>

    output {
        File vcf = "~{sample}.vcf"
    }
}

# pisces tumor only
# https://github.com/Illumina/Pisces
# 建立索引
# dotnet CreateGenomeSizeFile.dll \
#     -g hg19/ \
#     -s "Homo sapiens (UCSC hg19)" \
#     -o hg19/
task Pisces {
    input {
        String sample
        File bam
        File bai
        File bed
        Int threads
        Int depth
        Float maf
    }

    File PISCES = "/home/bioinfo/ubuntu/software/Pisces_5.2.10.49/Pisces.dll"
    File DOTNET = "/usr/bin/dotnet"
    File reference = "/home/bioinfo/ubuntu/database/hg19/ucsc.hg19.fasta"
    File ref_xml = "/home/bioinfo/ubuntu/database/hg19/GenomeSize.xml"
    String database = "/home/bioinfo/ubuntu/database/hg19/"

    command <<<
        ~{DOTNET} ~{PISCES} -b ~{bam} \
            -g ~{database} \
            -o ~{sample} \
            -t ~{threads} \
            -i ~{bed} \
            --mindp ~{depth} \
            --minvf ~{maf} \
            --minvq 0 --threadbychr true
    >>>

    output {
        File vcf = "~{sample}/{sample}.genome.vcf"
    }
}

task PiscesFix {
    input {
        String sample
        File vcf
    }

    File BCFTOOLS = "/home/bioinfo/ubuntu/software/bcftools-1.11/bcftools"

    command <<<
        ~{BCFTOOLS} view \
            -e "GT='0/0 | GT='./.' | GT='0/.'" \
            ~{vcf} > ~{sample}.muts.vcf
        ~{BCFTOOLS} view \
            -e "FILTER='LowDP'" \
            ~{sample}.muts.vcf > ~{sample}.pisces.vcf
    >>>

    output {
        File filterVcf = "~{sample}.pisces.vcf"
    }
}

# freebayes
# https://github.com/freebayes/freebayes
# 为了避免后续annovar注释的bug，此处设定genotyping-max-banddepth
# 同一位点只输出最多突变条数的突变方向
task Freebayes {
    input {
        String sample
        File bam
        File bai
        File bed
        Int depth
        Float maf
    }

    File reference = "/home/bioinfo/ubuntu/database/hg19/ucsc.hg19.fasta"
    File FREEBAYES = "/home/bioinfo/ubuntu/software/freebayes/freebayes"

    command <<<
        ~{FREEBAYES} -f ~{reference} \
            ~{bam} \
            -t ~{bed} \
            -F ~{maf} \
            -C 5 \
            --min-coverage ~{depth} \
            --genotyping-max-banddepth 1 \
            > ~{sample}.freebayes.vcf
        sed -i 's/0\/0/0\/1/g' ~{sample}.freebayes.vcf
    >>>

    output {
        File vcf = "~{sample}.freebayes.vcf"
    }
}


# 过滤
# bcftools
# http://samtools.github.io/bcftools/bcftools.html
task Filter {
    input {
        String sample
        File vcf
        Int depth
    }

    File BCFTOOLS = "/home/bioinfo/ubuntu/software/bcftools-1.11/bcftools"

    command <<<
        ~{BCFTOOLS} view \
            -i 'MIN(FORMAT/DP)>=~{depth}' \
            ~{vcf} \
            > ~{sample}.filter.vcf
    >>>

    output {
        File filterVcf = "~{sample}.filter.vcf"
    }
}

task FilterPair {
    input {
        String sample
        File vcf
        Int depth
    }

    File BCFTOOLS = "/home/bioinfo/ubuntu/software/bcftools-1.11/bcftools"

    command <<<
        ~{BCFTOOLS} view \
            -i 'FMT/DP[0]>=~{depth} && FMT/DP[1]>=100' \
            -s ~{sample} \
            ~{vcf} \
            > ~{sample}.filter.vcf
    >>>

    output {
        File filterVcf = "~{sample}.filter.vcf"
    }
}
