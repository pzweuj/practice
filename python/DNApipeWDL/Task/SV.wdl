version 1.0
# pzw
# 20210518
# SV


task Lumpy {
    input {
        String sample
        Int threads
        File bam
        File bai
    }

    File LUMPYEXPRESS = "/home/bioinfo/ubuntu/software/lumpy-sv/scripts/lumpyexpress"
    File EXTRACTSR = "/home/bioinfo/ubuntu/software/lumpy-sv/scripts/extractSplitReads_BwaMem"
    File SAMTOOLS = "/home/bioinfo/ubuntu/software/samtools-1.11/samtools"

    command <<<
        ~{SAMTOOLS} view -bh -F 1294 ~{bam} \
            | ~{SAMTOOLS} sort -@ ~{threads} - \
            -o ~{sample}.discordants.bam
        ~{SAMTOOLS} index ~{sample}.discordants.bam
        ~{SAMTOOLS} view -h ~{bam} \
            | ~{EXTRACTSR} \
            -i stdin \
            | ~{SAMTOOLS} view -bSh - \
            | ~{SAMTOOLS} sort -@ ~{threads} - \
            -o ~{sample}.splitters.bam
        ~{SAMTOOLS} index ~{sample}.splitters.bam
        ~{LUMPYEXPRESS} -B ~{bam} \
            -D ~{sample}.discordants.bam \
            -S ~{sample}.splitters.bam \
            -o ~{sample}.lumpy.vcf
    >>>

    output {
        File vcf = "~{sample}.lumpy.vcf"
    }
}


task LumpyPair {
    input {
        String sample
        String pair
        Int threads
        File tumorBam
        File tumorBamBai
        File normalBam
        File normalBamBai
    }

    File LUMPYEXPRESS = "/home/bioinfo/ubuntu/software/lumpy-sv/scripts/lumpyexpress"
    File EXTRACTSR = "/home/bioinfo/ubuntu/software/lumpy-sv/scripts/extractSplitReads_BwaMem"
    File SAMTOOLS = "/home/bioinfo/ubuntu/software/samtools-1.11/samtools"

    command <<<
        ~{SAMTOOLS} view -bh -F 1294 ~{tumorBam} \
            | ~{SAMTOOLS} sort -@ ~{threads} - \
            -o ~{sample}.discordants.bam
        ~{SAMTOOLS} index ~{sample}.discordants.bam
        ~{SAMTOOLS} view -h ~{tumorBam} \
            | ~{EXTRACTSR} \
            -i stdin \
            | ~{SAMTOOLS} view -bSh - \
            | ~{SAMTOOLS} sort -@ ~{threads} - \
            -o ~{sample}.splitters.bam
        ~{SAMTOOLS} index ~{sample}.splitters.bam
        ~{SAMTOOLS} view -bh -F 1294 ~{normalBam} \
            | ~{SAMTOOLS} sort -@ ~{threads} - \
            -o ~{pair}.discordants.bam
        ~{SAMTOOLS} index ~{pair}.discordants.bam
        ~{SAMTOOLS} view -h ~{normalBam} \
            | ~{EXTRACTSR} \
            -i stdin \
            | ~{SAMTOOLS} view -bSh - \
            | ~{SAMTOOLS} sort -@ ~{threads} - \
            -o ~{pair}.splitters.bam
        ~{SAMTOOLS} index ~{pair}.splitters.bam
        ~{LUMPYEXPRESS} -B ~{tumorBam},~{normalBam} \
            -D ~{sample}.discordants.bam,~{pair}.discordants.bam \
            -S ~{sample}.splitters.bam,~{pair}.splitters.bam \
            -o ~{sample}.lumpy.vcf
    >>>

    output {
        File vcf = "~{sample}.lumpy.vcf"
    }
}


task LumpyFilter {
    input {
        String sample
        File bam
        File bai
        File vcf
        Int depth
        Float maf
    }

    File SAMTOOLS = "/home/bioinfo/ubuntu/software/samtools-1.11/samtools"

    command <<<
        python3 <<CODE
        import os
        sampleID = "~{sample}"
        bamFile = "~{bam}"
        lumpyResultsFile = "~{vcf}"

        filterDP = ~{depth}
        filterVAF = ~{maf}
        filterAlt = int(filterDP * filterVAF)

        lumpyResults = open(lumpyResultsFile, "r")
        lumpyResultsFilter = open(lumpyResultsFile.replace(".vcf", ".filter.vcf"), "w")
        for line in lumpyResults:
            if line.startswith("#"):
                lumpyResultsFilter.write(line)
            else:
                if not "SVTYPE=BND" in line:
                    continue

                lineAfterSplit = line.split("\n")[0].split("\t")
                chrom1 = lineAfterSplit[0]
                pos1 = lineAfterSplit[1]
                alt = lineAfterSplit[4]
                if alt[0:2] == "N[":
                    chrom2s = alt.replace("N[", "").replace("[", "")
                elif alt[0:2] == "N]":
                    chrom2s = alt.replace("N]", "").replace("]", "")
                elif alt[-2:] == "[N":
                    chrom2s = alt.replace("[N", "").replace("[", "")
                elif alt[-2:] == "]N":
                    chrom2s = alt.replace("]N", "").replace("]", "")
                else:
                    continue

                chrom2 = chrom2s.split(":")[0]
                pos2 = chrom2s.split(":")[1]

                FORMAT = "PR:SR"
                FORMAT_origin = lineAfterSplit[8].split(":")
                FORMAT_info = lineAfterSplit[9]
                infos = FORMAT_info.split(":")
                formatDict = {}
                for F in range(len(FORMAT_origin)):
                    formatDict[FORMAT_origin[F]] = infos[F]

                try:
                    PR_alt = formatDict["PE"]
                except:
                    PR_alt = "0"
                try:
                    SR_alt = formatDict["SR"]
                except:
                    SR_alt = "0"
                check1 = chrom1 + ":" + pos1 + "-" + pos1
                check2 = chrom2 + ":" + pos2 + "-" + pos2
                
                # 突变reads数过滤
                if (int(PR_alt) + int(SR_alt)) < filterAlt:
                    continue

                dp1 = os.popen("~{SAMTOOLS} depth {bamFile} -r {check1}".format(bamFile=bamFile, check1=check1))
                dp2 = os.popen("~{SAMTOOLS} depth {bamFile} -r {check2}".format(bamFile=bamFile, check2=check2))
                for d1 in dp1.readlines():
                    if d1.startswith("chr"):
                        depth1 = d1.replace("\n", "").split("\t")[2]
                for d2 in dp2.readlines():
                    if d2.startswith("chr"):
                        depth2 = d2.replace("\n", "").split("\t")[2]
                DP = int(depth1) + int(depth2)

                # 深度过滤
                if DP < filterDP:
                    continue

                print(line)
                PR_ref = str(DP)
                SR_ref = "0"
                output = "\t".join(lineAfterSplit[0:8]) + "\t" + FORMAT + "\t" + PR_ref + "," + PR_alt + ":" + SR_ref + "," + SR_alt + "\n"
                lumpyResultsFilter.write(output)
        lumpyResultsFilter.close()
        lumpyResults.close()
        CODE
    >>>

    output {
        File filterVcf = "~{sample}.lumpy.filter.vcf"
    }

}

task LumpyFilterPair {
    input {
        String sample
        String pair
        File tumorBam
        File tumorBai
        File normalBam
        File normalBai
        File vcf
        Int depth
        Float maf
    }

    File SAMTOOLS = "/home/bioinfo/ubuntu/software/samtools-1.11/samtools"

    command <<<
        python3 <<CODE
        import os
        sampleID = "~{sample}"
        pairID = "~{pair}"
        bamFile = "~{tumorBam}"
        pairFile = "~{normalBam}"

        lumpyResultsFile = "~{vcf}"

        filterDP = ~{depth}
        filterVAF = ~{maf}
        filterAlt = int(filterDP * filterVAF)

        lumpyResults = open(lumpyResultsFile, "r")
        lumpyResultsFilter = open(lumpyResultsFile.replace(".vcf", ".filter.vcf"), "w")
        for line in lumpyResults:
            # print(line)
            if line.startswith("#"):
                line.replace(sampleID + "\t" + pairID, pairID + "\t" + sampleID)
                lumpyResultsFilter.write(line)
            else:
                if not "SVTYPE=BND" in line:
                    continue

                lineAfterSplit = line.split("\n")[0].split("\t")
                chrom1 = lineAfterSplit[0]
                pos1 = lineAfterSplit[1]
                alt = lineAfterSplit[4]
                if alt[0:2] == "N[":
                    chrom2s = alt.replace("N[", "").replace("[", "")
                elif alt[0:2] == "N]":
                    chrom2s = alt.replace("N]", "").replace("]", "")
                elif alt[-2:] == "[N":
                    chrom2s = alt.replace("[N", "").replace("[", "")
                elif alt[-2:] == "]N":
                    chrom2s = alt.replace("]N", "").replace("]", "")
                else:
                    continue

                chrom2 = chrom2s.split(":")[0]
                pos2 = chrom2s.split(":")[1]

                FORMAT = "PR:SR"
                FORMAT_origin = lineAfterSplit[8].split(":")
                FORMAT_info = lineAfterSplit[9]
                infos = FORMAT_info.split(":")
                formatDict = {}
                for F in range(len(FORMAT_origin)):
                    formatDict[FORMAT_origin[F]] = infos[F]
                
                if pairID != None:
                    FORMAT_info_n = lineAfterSplit[10]
                    infos_n = FORMAT_info_n.split(":")
                    for F in range(len(FORMAT_origin)):
                        formatDict_n[FORMAT_origin[F]] = infos_n[F]
                    try:
                        PR_alt_n = formatDict_n["PE"]
                    except:
                        PR_alt_n = "0"
                    try:
                        SR_alt_n = formatDict_n["SR"]
                    except:
                        SR_alt_n = "0"

                try:
                    PR_alt = formatDict["PE"]
                except:
                    PR_alt = "0"
                try:
                    SR_alt = formatDict["SR"]
                except:
                    SR_alt = "0"
                check1 = chrom1 + ":" + pos1 + "-" + pos1
                check2 = chrom2 + ":" + pos2 + "-" + pos2
                
                # 突变reads数过滤
                if (int(PR_alt) + int(SR_alt)) < filterAlt:
                    continue

                dp1 = os.popen("~{SAMTOOLS} depth {bamFile} -r {check1}".format(bamFile=bamFile, check1=check1))
                dp2 = os.popen("~{SAMTOOLS} depth {bamFile} -r {check2}".format(bamFile=bamFile, check2=check2))
                for d1 in dp1.readlines():
                    if d1.startswith("chr"):
                        depth1 = d1.replace("\n", "").split("\t")[2]
                for d2 in dp2.readlines():
                    if d2.startswith("chr"):
                        depth2 = d2.replace("\n", "").split("\t")[2]
                DP = int(depth1) + int(depth2)
                
                
                if pairID != None:
                    dp1_n = os.popen("~{SAMTOOLS} depth {pairFile} -r {check1}".format(pairFile=pairFile, check1=check1))
                    dp2_n = os.popen("~{SAMTOOLS} depth {pairFile} -r {check2}".format(pairFile=pairFile, check2=check2))
                    depth1_n = depth2_n = "0"
                    for d1_n in dp1_n.readlines():
                        if d1_n.startswith("chr"):
                            depth1_n = d1_n.replace("\n", "").split("\t")[2]
                    for d2_n in dp2_n.readlines():
                        if d2_n.startswith("chr"):
                            depth2_n = d2_n.replace("\n", "").split("\t")[2]
                    DP_n = int(depth1_n) + int(depth2_n)                    


                # 深度过滤
                if DP < filterDP:
                    continue

                print(line)
                PR_ref = str(DP)
                SR_ref = "0"

                if pairID != None:
                    PR_ref_n = str(DP_n)
                    SR_ref_n = "0"                    

                output = "\t".join(lineAfterSplit[0:8]) + "\t" + FORMAT + "\t" + PR_ref_n + "," + PR_alt_n + ":" + SR_ref_n + "," + SR_alt_n + "\t" + PR_ref + "," + PR_alt + ":" + SR_ref + "," + SR_alt + "\n"
                lumpyResultsFilter.write(output)
        lumpyResultsFilter.close()
        lumpyResults.close()
        CODE
    >>>

    output {
        File vcfFilter = "~{sample}.lumpy.filter.vcf"
    }

}

# manta
# https://github.com/Illumina/manta
# 当前建议使用manta，因为manta可同时输出read depth
task Manta {
    input {
        String sample
        Int threads
        File bam
        File bai
    }

    File reference = "/home/bioinfo/ubuntu/database/hg19/ucsc.hg19.fasta"
    File ref_fai = "/home/bioinfo/ubuntu/database/hg19/ucsc.hg19.fasta.fai"
    File MANTA = "/home/bioinfo/ubuntu/software/manta-1.6.0/bin/configManta.py"

    command <<<
        mkdir ~{sample}_tmp
        python3 ~{MANTA} \
            --tumorBam ~{bam} \
            --referenceFasta ~{reference} \
            --exome \
            --generateEvidenceBam \
            --runDir ~{sample}_tmp
        ~{sample}_tmp/runWorkflow.py -j ~{threads}
        zcat ~{sample}_tmp/results/variants/tumorSV.vcf.gz > ~{sample}.manta.vcf
    >>>

    output {
        File vcf = "~{sample}.manta.vcf"
    }
}

task MantaPair {
    input {
        String sample
        Int threads
        File tumorBam
        File tumorBai
        File normalBam
        File normalBai
    }

    File reference = "/home/bioinfo/ubuntu/database/hg19/ucsc.hg19.fasta"
    File ref_fai = "/home/bioinfo/ubuntu/database/hg19/ucsc.hg19.fasta.fai"
    File MANTA = "/home/bioinfo/ubuntu/software/manta-1.6.0/bin/configManta.py"

    command <<<
        mkdir ~{sample}_tmp
        python3 ~{MANTA} \
            --tumorBam ~{tumorBam} \
            --normalBam ~{normalBam} \
            --referenceFasta ~{reference} \
            --exome \
            --generateEvidenceBam \
            --runDir ~{sample}_tmp
        ~{sample}_tmp/runWorkflow.py -j ~{threads}
        zcat ~{sample}_tmp/results/variants/somaticSV.vcf.gz > ~{sample}.manta.vcf        
    >>>

    output {
        File vcf = "~{sample}.manta.vcf"
    }

}

task MantaFilter {
    input {
        String sample
        File vcf
    }

    
}