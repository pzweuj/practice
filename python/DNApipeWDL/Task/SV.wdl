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

    command <<<
        samtools view -bh -F 1294 ~{bam} \
            | samtools sort -@ ~{threads} - \
            -o ~{sample}.discordants.bam
        samtools index ~{sample}.discordants.bam
        samtools view -h ~{bam} \
            | extractSplitReads_BwaMem \
            -i stdin \
            | samtools view -bSh - \
            | samtools sort -@ ~{threads} - \
            -o ~{sample}.splitters.bam
        samtools index ~{sample}.splitters.bam
        lumpyexpress -B ~{bam} \
            -D ~{sample}.discordants.bam \
            -S ~{sample}.splitters.bam \
            -o ~{sample}.lumpy.vcf
    >>>

    output {
        File vcf = "~{sample}.lumpy.vcf"
    }

    runtime {
        docker: "marrip/lumpy:v0.3.1"
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

    command <<<
        samtools view -bh -F 1294 ~{tumorBam} \
            | samtools sort -@ ~{threads} - \
            -o ~{sample}.discordants.bam
        samtools index ~{sample}.discordants.bam
        samtools view -h ~{tumorBam} \
            | extractSplitReads_BwaMem \
            -i stdin \
            | samtools view -bSh - \
            | samtools sort -@ ~{threads} - \
            -o ~{sample}.splitters.bam
        samtools} index ~{sample}.splitters.bam
        samtools view -bh -F 1294 ~{normalBam} \
            | samtools sort -@ ~{threads} - \
            -o ~{pair}.discordants.bam
        samtools index ~{pair}.discordants.bam
        samtools view -h ~{normalBam} \
            | extractSplitReads_BwaMem \
            -i stdin \
            | samtools view -bSh - \
            | samtools sort -@ ~{threads} - \
            -o ~{pair}.splitters.bam
        samtools index ~{pair}.splitters.bam
        lumpyexpress -B ~{tumorBam},~{normalBam} \
            -D ~{sample}.discordants.bam,~{pair}.discordants.bam \
            -S ~{sample}.splitters.bam,~{pair}.splitters.bam \
            -o ~{sample}.lumpy.vcf
    >>>

    output {
        File vcf = "~{sample}.lumpy.vcf"
    }

    runtime {
        docker: "marrip/lumpy:v0.3.1"
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
        lumpyResultsFilter = open("~{sample}.lumpy.filter.vcf", "w")
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
        File vcfFilter = "~{sample}.lumpy.filter.vcf"
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
        lumpyResultsFilter = open("~{sample}.lumpy.filter.vcf", "w")
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

    command <<<
        mkdir ~{sample}_tmp
        configManta.py \
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

    runtime {
        docker: "pzweuj/manta:v1.6.0"
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

    command <<<
        mkdir ~{sample}_tmp
        configManta.py \
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

    runtime {
        docker: "pzweuj/manta:v1.6.0"
    }
}

task MantaFilter {
    input {
        String sample
        File vcf
        Int depth
        Float maf
    }

    command <<<
        python3 <<CODE
        sampleID = "~{sample}"
        mantaResultsFile = "~{vcf}"

        filterDP = ~{depth}
        filterVAF = ~{maf}
        filterAlt = int(filterDP * filterVAF)

        mantaResults = open(mantaResultsFile, "r")
        mantaResultsFilter = open("~{sample}.manta.filter.vcf", "w")
        for line in mantaResults:
            if line.startswith("#"):
                mantaResultsFilter.write(line)
            else:
                lineAfterSplit = line.split("\n")[0].split("\t")
                FORMAT = lineAfterSplit[8]
                
                FORMAT_info = lineAfterSplit[9]
                if "PR" not in FORMAT:
                    SR = FORMAT_info.split(",")
                    PR_ref = "0"
                    PR_alt = "0"
                    SR_ref = SR[0]
                    SR_alt = SR[1]
                elif "SR" not in FORMAT:
                    PR = FORMAT_info.split(",")
                    SR_ref = "0"
                    SR_alt = "0"
                    PR_ref = PR[0]
                    PR_alt = PR[1]
                else:
                    FORMAT_infos = FORMAT_info.split(":")
                    PR = FORMAT_infos[0].split(",")
                    SR = FORMAT_infos[1].split(",")
                    PR_ref = PR[0]
                    PR_alt = PR[1]
                    SR_ref = SR[0]
                    SR_alt = SR[1]
                
                Ref_total = int(PR_ref) + int(SR_ref)
                Alt_total = int(PR_alt) + int(SR_alt)
                if Ref_total < filterDP:
                    continue
                if Alt_total < filterAlt:
                    continue

                mantaResultsFilter.write(line)
        mantaResultsFilter.close()
        mantaResults.close()        
        CODE
    >>>

    output {
        File vcfFilter = "~{sample}.manta.filter.vcf"
    }
}


task MantaFilterPair {
    input {
        String sample
        String pair
        File vcf
        Int depth
        Float maf
    }

    command <<<
        sampleID = "~{sample}"
        pairID = "~{pair}"
        mantaResultsFile = "~{vcf}"

        filterDP = ~{depth}
        filterVAF = ~{maf}
        filterAlt = int(filterDP * filterVAF)

        mantaResults = open(mantaResultsFile, "r")
        mantaResultsFilter = open("~{sample}.manta.filter.vcf", "w")
        for line in mantaResults:
            if line.startswith("#"):
                mantaResultsFilter.write(line)
            else:
                lineAfterSplit = line.split("\n")[0].split("\t")
                FORMAT = lineAfterSplit[8]
                FORMAT_info = lineAfterSplit[10]

                if "PR" not in FORMAT:
                    SR = FORMAT_info.split(",")
                    PR_ref = "0"
                    PR_alt = "0"
                    SR_ref = SR[0]
                    SR_alt = SR[1]
                elif "SR" not in FORMAT:
                    PR = FORMAT_info.split(",")
                    SR_ref = "0"
                    SR_alt = "0"
                    PR_ref = PR[0]
                    PR_alt = PR[1]
                else:
                    FORMAT_infos = FORMAT_info.split(":")
                    PR = FORMAT_infos[0].split(",")
                    SR = FORMAT_infos[1].split(",")
                    PR_ref = PR[0]
                    PR_alt = PR[1]
                    SR_ref = SR[0]
                    SR_alt = SR[1]
                
                Ref_total = int(PR_ref) + int(SR_ref)
                Alt_total = int(PR_alt) + int(SR_alt)
                if Ref_total < filterDP:
                    continue
                if Alt_total < filterAlt:
                    continue

                mantaResultsFilter.write(line)
        mantaResultsFilter.close()
        mantaResults.close()   
    >>>

    output {
        File vcfFilter = "~{sample}.manta.filter.vcf"
    }
}


task FusionAnno {
    input {
        String sample
        String pair = ""
        File lumpyVcf
        File mantaVcf
    }

    File refFlat = "/home/bioinfo/ubuntu/database/hg19/hg19_refFlat.txt"

    command <<<
        python3 <<CODE
        import os
        def checkFusionHotSpot(fusionStrand, chrom, breakPoint):
            database = "~{refFlat}"
            db = open(database, "r")
            output = []
            ts_output = []
            for line in db:
                l = line.split("\n")[0].split("\t")
                gene = l[0]
                ts = l[1]
                chrom_db = l[2]
                strand = l[3]
                gene_start = l[4]
                gene_end = l[5]
                cds_start = l[6]
                cds_end = l[7]
                exon_nums = l[8]
                exon_starts = l[9]
                exon_ends = l[10]

                # 寻找注释
                if chrom == chrom_db:
                    # 找到目标基因
                    if (int(breakPoint) >= int(gene_start)) and (int(breakPoint) <= int(gene_end)):

                        exon_starts_list = exon_starts.split(",")[0: -1]
                        exon_ends_list = exon_ends.split(",")[0: -1]

                        # 将外显子起止点形成列表
                        exon_checkpoint = []
                        i = 0
                        while i < int(exon_nums):
                            exon_checkpoint.append(int(exon_starts_list[i]))
                            exon_checkpoint.append(int(exon_ends_list[i]))
                            i += 1

                        # 判断目标基因的转录方向
                        if strand == "+":
                            if int(breakPoint) < exon_checkpoint[0]:
                                if fusionStrand == "+":
                                    exon_out = gene + "_5UTR"
                                elif fusionStrand == "-":
                                    exon_out = gene + "_exon1"
                                else:
                                    exon_out = gene + "_?"

                            elif int(breakPoint) > exon_checkpoint[-1]:
                                if fusionStrand == "+":
                                    exon_out = gene + "_exon" + exon_nums
                                elif fusionStrand == "-":
                                    exon_out = gene + "_3UTR"
                                else:
                                    exon_out = gene + "_?"

                            else:
                                n = 0
                                while n < (2 * int(exon_nums) - 1):
                                    if (int(breakPoint) >= exon_checkpoint[n]) and (int(breakPoint) <= exon_checkpoint[n+1]):
                                        checkN = n
                                    n += 1

                                if fusionStrand == "+":
                                    exon_out = gene + "_exon" + str(checkN // 2 + 1)
                                elif fusionStrand == "-":
                                    exon_out = gene + "_exon" + str((checkN + 1) // 2 + 1)
                                else:
                                    exon_out = gene + "_?"


                        if strand == "-":
                            if int(breakPoint) < exon_checkpoint[0]:
                                if fusionStrand == "+":
                                    exon_out = gene + "_3UTR"
                                elif fusionStrand == "-":
                                    exon_out = gene + "_exon" + exon_nums
                                else:
                                    exon_out = gene + "_?"


                            elif int(breakPoint) > exon_checkpoint[-1]:
                                if fusionStrand == "+":
                                    exon_out = gene + "_exon1"
                                elif fusionStrand == "-":
                                    exon_out = gene + "_5UTR"
                                else:
                                    exon_out = gene + "_?"

                            else:
                                n = 0
                                while n < (2 * int(exon_nums) - 1):
                                    if (int(breakPoint) >= exon_checkpoint[n]) and (int(breakPoint) <= exon_checkpoint[n+1]):
                                        checkN = n
                                    n += 1

                                if fusionStrand == "+":
                                    exon_out = gene + "_exon" + str(int(exon_nums) - checkN // 2)
                                elif fusionStrand == "-":
                                    exon_out = gene + "_exon" + str(int(exon_nums) - (checkN + 1) // 2)
                                else:
                                    exon_out = gene + "_?"

                        output.append(exon_out)
                        ts_output.append(ts)
            if not output:
                output.append("GeneUnknown_exon?")        
            if not ts_output:
                ts_output.append("Transcript?")
            db.close()
            return [output, ts_output]

        # 注释开始
        sampleID = "~{sample}"
        pairID = "~{pair}"

        # lumpy
        svFile = open("~{lumpyVcf}", "r")
        svAnno = open("~{sample}.fusion.txt", "w")
        svAnno.write("chrom1\tbreakpoint1\tgene1\tchrom2\tbreakpoint2\tgene2\tfusionType\tAlt\tgeneSymbol\tPR\tSR\tDP\tVAF\tExon\tTranscript\tSource\n")
        for line in svFile:
            if line.startswith("#"):
                continue
            else:
                lineAfterSplit = line.split("\t")
                chrom = lineAfterSplit[0]
                Pos = lineAfterSplit[1]
                ID = lineAfterSplit[2]
                Ref = lineAfterSplit[3]
                Alt = lineAfterSplit[4]
                Qual = lineAfterSplit[5]
                Filter = lineAfterSplit[6]
                Info = lineAfterSplit[7]
                Format = lineAfterSplit[8]

                if pairID == "":
                    Sample = lineAfterSplit[9]
                else:
                    Sample = lineAfterSplit[10]

                if chrom == "chrM":
                    continue

                if ":" in Format:
                    S = Sample.split(":")
                    PR = S[0].split(",")
                    SR = S[1].split(",")
                    PR_ref = int(PR[0])
                    PR_alt = int(PR[1])
                    SR_ref = int(SR[0])
                    SR_alt = int(SR[1])
                else:
                    PR = Sample.split(",")
                    PR_ref = int(PR[0])
                    PR_alt = int(PR[1])
                    SR_ref = SR_alt = 0

                DP = PR_ref + SR_ref + PR_alt + SR_alt
                VAF = "%.2f" % (100 * float(PR_alt + SR_alt) / DP) + "%"

                if ("[" in Alt) or ("]" in Alt):
                    if chrom in Alt:
                        fusionType = "Inversion"
                    else:
                        fusionType = "Translocation"

                    # G    [chr1:1111111[G
                    if Alt[0] == "[":
                        chrom2 = Alt.split(":")[0].split("[")[1]
                        breakpoint2 = Alt.split(":")[1].split("[")[0]
                        strand = "--"
                    
                    # G    G[chr1:1111111[
                    elif Alt[-1] == "[":
                        chrom2 = Alt.split(":")[0].split("[")[1]
                        breakpoint2 = Alt.split(":")[1].split("[")[0]
                        strand = "+-"

                    # G    ]chr1:1111111]G
                    elif Alt[0] == "]":
                        chrom2 = Alt.split(":")[0].split("]")[1]
                        breakpoint2 = Alt.split(":")[1].split("]")[0]
                        strand = "-+"

                    # G    G]chr1:1111111]
                    elif Alt[-1] == "]":
                        chrom2 = Alt.split(":")[0].split("]")[1]
                        breakpoint2 = Alt.split(":")[1].split("]")[0]
                        strand = "++"

                    else:
                        chrom2 = "Unknown"
                        breakpoint2 = "Unknown"
                        strand = "Unknown"
                        print("can not analysis: " + ID)

                    anno1 = checkFusionHotSpot(strand[0], chrom, Pos)
                    if chrom2 != "Unknown":
                        anno2 = checkFusionHotSpot(strand[1], chrom2, breakpoint2)
                    else:
                        anno2 = [["GeneUnknown_exon?"], ["Transcript?"]]

                    gene1 = anno1[0][0].split("_")[0]
                    gene2 = anno2[0][0].split("_")[0]

                    if gene1 == gene2:
                        continue

                    exon_output = ",".join(anno1[0]) + "|" + ",".join(anno2[0])
                    transcript_output = ",".join(anno1[1]) + "|" + ",".join(anno2[1])
                    outputStringList = [chrom, Pos, gene1, chrom2, breakpoint2, gene2, fusionType, Alt, gene1 + "-" + gene2, str(PR_alt), str(SR_alt), str(DP), VAF, exon_output, transcript_output]

                    outputString = "\t".join(outputStringList)
                    print(outputString)
                    svAnno.write(outputString + "\tLumpy\n")
        svFile.close()
        
        # manta
        svFile = open("~{mantaVcf}", "r")
        for line in svFile:
            if line.startswith("#"):
                continue
            else:
                lineAfterSplit = line.split("\t")
                chrom = lineAfterSplit[0]
                Pos = lineAfterSplit[1]
                ID = lineAfterSplit[2]
                Ref = lineAfterSplit[3]
                Alt = lineAfterSplit[4]
                Qual = lineAfterSplit[5]
                Filter = lineAfterSplit[6]
                Info = lineAfterSplit[7]
                Format = lineAfterSplit[8]

                if pairID == "":
                    Sample = lineAfterSplit[9]
                else:
                    Sample = lineAfterSplit[10]

                if chrom == "chrM":
                    continue

                if ":" in Format:
                    S = Sample.split(":")
                    PR = S[0].split(",")
                    SR = S[1].split(",")
                    PR_ref = int(PR[0])
                    PR_alt = int(PR[1])
                    SR_ref = int(SR[0])
                    SR_alt = int(SR[1])
                else:
                    PR = Sample.split(",")
                    PR_ref = int(PR[0])
                    PR_alt = int(PR[1])
                    SR_ref = SR_alt = 0

                DP = PR_ref + SR_ref + PR_alt + SR_alt
                VAF = "%.2f" % (100 * float(PR_alt + SR_alt) / DP) + "%"

                if ("[" in Alt) or ("]" in Alt):
                    if chrom in Alt:
                        fusionType = "Inversion"
                    else:
                        fusionType = "Translocation"

                    # G    [chr1:1111111[G
                    if Alt[0] == "[":
                        chrom2 = Alt.split(":")[0].split("[")[1]
                        breakpoint2 = Alt.split(":")[1].split("[")[0]
                        strand = "--"
                    
                    # G    G[chr1:1111111[
                    elif Alt[-1] == "[":
                        chrom2 = Alt.split(":")[0].split("[")[1]
                        breakpoint2 = Alt.split(":")[1].split("[")[0]
                        strand = "+-"

                    # G    ]chr1:1111111]G
                    elif Alt[0] == "]":
                        chrom2 = Alt.split(":")[0].split("]")[1]
                        breakpoint2 = Alt.split(":")[1].split("]")[0]
                        strand = "-+"

                    # G    G]chr1:1111111]
                    elif Alt[-1] == "]":
                        chrom2 = Alt.split(":")[0].split("]")[1]
                        breakpoint2 = Alt.split(":")[1].split("]")[0]
                        strand = "++"

                    else:
                        chrom2 = "Unknown"
                        breakpoint2 = "Unknown"
                        strand = "Unknown"
                        print("can not analysis: " + ID)

                    anno1 = checkFusionHotSpot(strand[0], chrom, Pos)
                    if chrom2 != "Unknown":
                        anno2 = checkFusionHotSpot(strand[1], chrom2, breakpoint2)
                    else:
                        anno2 = [["GeneUnknown_exon?"], ["Transcript?"]]

                    gene1 = anno1[0][0].split("_")[0]
                    gene2 = anno2[0][0].split("_")[0]

                    if gene1 == gene2:
                        continue

                    exon_output = ",".join(anno1[0]) + "|" + ",".join(anno2[0])
                    transcript_output = ",".join(anno1[1]) + "|" + ",".join(anno2[1])
                    outputStringList = [chrom, Pos, gene1, chrom2, breakpoint2, gene2, fusionType, Alt, gene1 + "-" + gene2, str(PR_alt), str(SR_alt), str(DP), VAF, exon_output, transcript_output]

                    outputString = "\t".join(outputStringList)
                    print(outputString)
                    svAnno.write(outputString + "\tManta\n")        
        svFile.close()
        svAnno.close()

        CODE
    >>>

    output {
        File svAnno = "~{sample}.fusion.txt"
    }
}