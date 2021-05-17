version 1.0
# pzweuj
# 20210517
# QC

# fastp
# https://github.com/OpenGene/fastp
task Fastp {
    input {
        String sample
        File rawRead1
        File rawRead2
        Int threads        
    }

    File FASTP = "/home/bioinfo/ubuntu/software/fastp/fastp"

    command <<<
        ~{FASTP} \
            -i ~{rawRead1} \
            -I ~{rawRead2} \
            -o ~{sample}.clean_R1.fastq.gz \
            -O ~{sample}.clean_R2.fastq.gz \
            -w ~{threads} \
            -j ~{sample}.json \
            -h ~{sample}.html
    >>>

    output {
        File cleanRead1 = "~{sample}.clean_R1.fastq.gz"
        File cleanRead2 = "~{sample}.clean_R2.fastq.gz"
        File jsonReport = "~{sample}.json"
        File htmlReport = "~{sample}.html"
    }
}


# fastp umi模式
# https://github.com/OpenGene/fastp
task FastpUmi {
    input {
        String sample
        File rawRead1
        File rawRead2
        Int threads
        String umiLoc
        Int umiLen
    }

    File FASTP = "/home/bioinfo/ubuntu/software/fastp/fastp"

    command <<<
        ~{FASTP} \
            -i ~{rawRead1} \
            -I ~{rawRead2} \
            -o ~{sample}.clean_R1.fastq.gz \
            -O ~{sample}.clean_R2.fastq.gz \
            -w ~{threads} \
            -j ~{sample}.json \
            -h ~{sample}.html \
            -A -U --umi_prefix=UMI --umi_loc=~{umiLoc} --umi_len=~{umiLen} --umi_skip=2
    >>>

    output {
        File cleanRead1 = "~{sample}.clean_R1.fastq.gz"
        File cleanRead2 = "~{sample}.clean_R2.fastq.gz"
        File jsonReport = "~{sample}.json"
        File htmlReport = "~{sample}.html"
    }
}

# fastp结果整理
task FastpFormat {
    input {
        File jsonReport
        String sample
    }

    command <<<
        python3 <<CODE
        import json
        fastpJsonReport = "~{jsonReport}"
        jsonFile = json.load(open(fastpJsonReport, "r"))
        rawReads = jsonFile["summary"]["before_filtering"]["total_reads"]
        rawBases = jsonFile["summary"]["before_filtering"]["total_bases"]
        rawQ20Bases = jsonFile["summary"]["before_filtering"]["q20_bases"]
        rawQ30Bases = jsonFile["summary"]["before_filtering"]["q30_bases"]
        rawQ20Rate = "%.2f" % (float(rawQ20Bases / rawBases) * 100) + "%"
        rawQ30Rate = "%.2f" % (float(rawQ30Bases / rawBases) * 100) + "%"
        rawGC = "%.2f" % (jsonFile["summary"]["before_filtering"]["gc_content"] * 100) + "%"
        duplicationRate = "%.2f" % (jsonFile["duplication"]["rate"] * 100) + "%"
        cleanReads = jsonFile["summary"]["after_filtering"]["total_reads"]
        cleanBases = jsonFile["summary"]["after_filtering"]["total_bases"]
        cleanQ20Bases = jsonFile["summary"]["after_filtering"]["q20_bases"]
        cleanQ30Bases = jsonFile["summary"]["after_filtering"]["q30_bases"]
        cleanQ20Rate = "%.2f" % (float(cleanQ20Bases / cleanBases) * 100) + "%"
        cleanQ30Rate = "%.2f" % (float(cleanQ30Bases / cleanBases) * 100) + "%"
        cleanGC = "%.2f" % (jsonFile["summary"]["after_filtering"]["gc_content"] * 100) + "%"
        fastpReport = open("~{sample}.fastp.txt", "w")
        fastpReport.write("~{sample}" + " fastp QC Report\n")
        fastpReport.write("rawReads\t" + str(rawReads) + "\n")
        fastpReport.write("rawBases\t" + str(rawBases) + "\n")
        fastpReport.write("raw_q20\t" + str(rawQ20Bases) + "\n")
        fastpReport.write("raw_q20_rate\t" + rawQ20Rate + "\n")
        fastpReport.write("raw_q30\t" + str(rawQ30Bases) + "\n")
        fastpReport.write("raw_q30_rate\t" + rawQ30Rate + "\n")
        fastpReport.write("raw_GC_content\t" + rawGC + "\n")
        fastpReport.write("duplicationRate\t" + duplicationRate + "\n")
        fastpReport.write("cleanReads\t" + str(cleanReads) + "\n")
        fastpReport.write("cleanBases\t" + str(cleanBases) + "\n")
        fastpReport.write("clean_q20\t" + str(cleanQ20Bases) + "\n")
        fastpReport.write("clean_q20_rate\t" + cleanQ20Rate + "\n")
        fastpReport.write("clean_q30\t" + str(cleanQ30Bases) + "\n")
        fastpReport.write("clean_q30_rate\t" + cleanQ30Rate + "\n")
        fastpReport.write("clean_GC_content\t" + cleanGC + "\n")
        fastpReport.close()
        CODE
    >>>
    
    output {
        File fastpReport = "~{sample}.fastp.txt"
    }
}