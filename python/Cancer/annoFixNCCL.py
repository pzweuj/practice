# -*- coding:utf-8 -*-
# 20210910

import sys

annvarResultsFile = sys.argv[1]
annvarFixFile = sys.argv[2]
refTs = sys.argv[3]

def refTranscript(ts):
    refT = open(ts, "r")
    refDict = {}
    for r in refT:
        if not r.startswith("#"):
            gene = r.split("\t")[0]
            ts = r.split("\t")[1]
            refDict[gene] = ts
    refT.close()
    return refDict

def AATranslate(seq):
    transDict = {
        "Ala": "A", "Arg": "R", "Asn": "N",
        "Asp": "D", "Cys": "C", "Gln": "Q",
        "Glu": "E", "Gly": "G", "His": "H",
        "Ile": "I", "Leu": "L", "Lys": "K",
        "Met": "M", "Phe": "F", "Pro": "P",
        "Ser": "S", "Thr": "T", "Trp": "W",
        "Tyr": "Y", "Val": "V"
    }
    for i in transDict.keys():
        if i in seq:
            seq = seq.replace(i, transDict[i])
    
    return seq

# annvarResultsFile = "~{annovarResults}"
# annvarFixFile = "~{sample}.Anno.txt"

annvarResults = open(annvarResultsFile, "r")
results = open(annvarFixFile, "w")
results.write("Chr\tStart\tEnd\tRef\tAlt\tGene\tType\tTranscript\tcHGVS\tpHGVS\tVAF\tConsequence\t"\
    "AffectedExon\tDepth\tAlt_AD\tFunc.refGene\tAAChange.refGene\tFilterInfo\tavsnp150\tAF\tAF_popmax\tAF_male\t"\
    "AF_female\tAF_eas\t1000g2015aug_all\tExAC_ALL\tCLNDN\tCLNSIG\tInterVar_automated\tInterVar_sig\tJax_Ckb_Variants_Summary\tJax_Ckb_Drug_Summary\t"\
    "Civic\tOncoKB\tcosmic92_coding\tCADD_phred\tSIFT_pred\tPolyphen2_HDIV_pred\tPolyphen2_HVAR_pred\t"\
    "LRT_pred\tMutationTaster_pred\tMutationAssessor_pred\tFATHMM_pred\tPROVEAN_pred\t"\
    "M-CAP_pred\tREVEL_score\tdbscSNV_ADA_SCORE\tdbscSNV_RF_SCORE\tSpliceAI\n")

# 生成字典keys
lineDict = {}
lineDictList = []
for line in annvarResults:
    if line.startswith("Chr\tStart"):
        lines = line.replace("\n", "").split("\t")
        for k in lines:
            lineDict[k] = "-"
            lineDictList.append(k)
annvarResults.close()

# 再次打开
annvarResults = open(annvarResultsFile, "r")
for line in annvarResults:
    if not line.startswith("Chr\tStart"):
        lines = line.replace("\n", "").split("\t")
        
        # 生成字典
        for i in range(len(lines)):
            lineDict[lineDictList[i]] = lines[i]
        
        # 开始处理
        # 处理Insertion状态下的End
        if lineDict["Ref"] == "-":
            lineDict["Type"] = "Insertion"
            lineDict["End"] = "-"
        elif lineDict["Alt"] == "-":
            lineDict["Type"] = "Deletion"
        elif len(lineDict["Ref"]) == 1 and len(lineDict["Alt"]) == 1:
            lineDict["Type"] = "SNV"
        else:
            lineDict["Type"] = "Complex"

        # filter信息
        lineDict["FilterInfo"] = lineDict["Otherinfo10"]

        # 处理基因，使用snpeff的注释信息
        otif11 = lineDict["Otherinfo11"].split(";")
        snpeff_anno = "ANN=?|N|N|N|N|transcript|N|N|N|N|N|N|N|N|N|"
        for o11 in otif11:
            if o11.startswith("ANN="):
                snpeff_anno = o11
        snpeff_annos = snpeff_anno.replace("ANN=", "").split(",")
        
        # 此处需使用固定转录本
        if len(snpeff_annos) > 1:
            snpeff_select = snpeff_annos[0]
            gene = snpeff_select.split("|")[3]
            refTrans = refTranscript(refTs)
            try:
                ts = refTrans[gene]
                for sa in snpeff_annos:
                    if ts.split(".")[0] in sa:
                        if sa.split("|")[6] != "":
                            snpeff_select = sa
            except:
                snpeff_select = snpeff_select
        else:
            snpeff_select = snpeff_annos[0]
        s = snpeff_select.split("|")
        for ss in range(len(s)):
            if s[ss] == "":
                s[ss] = "-"

        lineDict["Gene"] = s[3]
        lineDict["Transcript"] = s[6]
        lineDict["cHGVS"] = s[9]
        lineDict["pHGVS"] = s[10]
        lineDict["AffectedExon"] = s[8]

        # 处理Consequence
        if lineDict["Func.refGene"] == "exonic" and lineDict["Type"] == "Complex":
            if "_" in lineDict["pHGVS"]:
                lineDict["Consequence"] = "Complex_mutation"
            elif "stop_gained" in s[1]:
                lineDict["Consequence"] = "Nonsense_substitution"
            elif "missense" in s[1]:
                lineDict["Consequence"] = "Missense_substitution"
            elif "synonymous" in s[1]:
                lineDict["Consequence"] = "Synonymous_substitution"
            else:
                lineDict["Consequence"] = "Complex_mutation"
        elif "missense" in s[1]:
            lineDict["Consequence"] = "Missense_substitution"
        elif "splice" in s[1]:
            lineDict["Consequence"] = "Splice_Site_mutation"
        elif "synonymous" in s[1]:
            lineDict["Consequence"] = "Synonymous_substitution"
        elif "inframe_deletion" in s[1]:
            lineDict["Consequence"] = "Inframe_deletion"
        elif "inframe_insertion" in s[1]:
            lineDict["Consequence"] = "Inframe_insertion"
        elif "frameshift" in s[1]:
            if lineDict["Type"] == "Insertion":
                lineDict["Consequence"] = "Frameshift_insertion"
            elif lineDict["Type"] == "Deletion":
                lineDict["Consequence"] = "Frameshift_deletion"
            else:
                lineDict["Consequence"] = "Frameshift_variant"
        elif "stop_gained" in s[1]:
            lineDict["Consequence"] = "Nonsense_substitution"
            if lineDict["pHGVS"] == "-":
                lineDict["Consequence"] = "Truncation_mutation"
        elif "stop_lost" in s[1]:
            lineDict["Consequence"] = "Elongation_mutation"
        else:
            lineDict["Consequence"] = "Other"

        # 处理VAF
        otif13 = lineDict["Otherinfo13"].split(":")
        otif12 = lineDict["Otherinfo12"].split(":")        
        format_zip = list(zip(otif12, otif13))
        format_dict = {}
        for f in format_zip:
            format_dict[f[0]] = f[1]
        lineDict["DP"] = format_dict["DP"]
        lineDict["Alt_AD"] = format_dict["AD"]
        if not "," in lineDict["Alt_AD"]:
            lineDict["Alt_AD"] = lineDict["Alt_AD"]
        else:
            lineDict["Alt_AD"] = lineDict["Alt_AD"].split(",")[1]
        lineDict["Ref_AD"] = str(int(lineDict["DP"]) - int(lineDict["Alt_AD"]))

        try:
            AF = "%.2f" % ((float(lineDict["Alt_AD"]) / float(lineDict["DP"])) * 100) + "%"
        except Exception:
            AF = "-"
        lineDict["VAF"] = AF

        # 合并证据等级
        lineDict["InterVar_sig"] = "PVS1=" + lineDict["PVS1"] + ";" + \
            "PS1=" + lineDict["PS1"] + ";" + \
            "PS2=" + lineDict["PS2"] + ";" + \
            "PS3=" + lineDict["PS3"] + ";" + \
            "PS4=" + lineDict["PS4"] + ";" + \
            "PM1=" + lineDict["PM1"] + ";" + \
            "PM2=" + lineDict["PM2"] + ";" + \
            "PM3=" + lineDict["PM3"] + ";" + \
            "PM4=" + lineDict["PM4"] + ";" + \
            "PM5=" + lineDict["PM5"] + ";" + \
            "PM6=" + lineDict["PM6"] + ";" + \
            "PP1=" + lineDict["PP1"] + ";" + \
            "PP2=" + lineDict["PP2"] + ";" + \
            "PP3=" + lineDict["PP3"] + ";" + \
            "PP4=" + lineDict["PP4"] + ";" + \
            "PP5=" + lineDict["PP5"] + ";" + \
            "BA1=" + lineDict["BA1"] + ";" + \
            "BS1=" + lineDict["BS1"] + ";" + \
            "BS2=" + lineDict["BS2"] + ";" + \
            "BS3=" + lineDict["BS3"] + ";" + \
            "BS4=" + lineDict["BS4"] + ";" + \
            "BP1=" + lineDict["BP1"] + ";" + \
            "BP2=" + lineDict["BP2"] + ";" + \
            "BP3=" + lineDict["BP3"] + ";" + \
            "BP4=" + lineDict["BP4"] + ";" + \
            "BP5=" + lineDict["BP5"] + ";" + \
            "BP6=" + lineDict["BP6"] + ";" + \
            "BP7=" + lineDict["BP7"]

        # 整理输出结果
        output = [lineDict["Chr"], lineDict["Start"], lineDict["End"], lineDict["Ref"],
            lineDict["Alt"], lineDict["Gene"], lineDict["Type"], lineDict["Transcript"],
            lineDict["cHGVS"], lineDict["pHGVS"], lineDict["VAF"], lineDict["Consequence"],
            lineDict["AffectedExon"], lineDict["DP"], lineDict["Alt_AD"], lineDict["Func.refGene"],
            lineDict["AAChange.refGene"], lineDict["FilterInfo"], lineDict["avsnp150"], lineDict["AF"], lineDict["AF_popmax"], lineDict["AF_male"], 
            lineDict["AF_female"], lineDict["AF_eas"], lineDict["1000g2015aug_all"], lineDict["ExAC_ALL"], lineDict["CLNDN"], lineDict["CLNSIG"],
            lineDict["InterVar_automated"], lineDict["InterVar_sig"], lineDict["Jax_Ckb_Variants_Summary"],
            lineDict["Jax_Ckb_Drug_Summary"], lineDict["Civic"], lineDict["OncoKB"], lineDict["cosmic92_coding"],
            lineDict["CADD_phred"], lineDict["SIFT_pred"], lineDict["Polyphen2_HDIV_pred"],
            lineDict["Polyphen2_HVAR_pred"], lineDict["LRT_pred"], lineDict["MutationTaster_pred"],
            lineDict["MutationAssessor_pred"], lineDict["FATHMM_pred"], lineDict["PROVEAN_pred"],
            lineDict["M-CAP_pred"], lineDict["REVEL_score"], lineDict["dbscSNV_ADA_SCORE"], lineDict["dbscSNV_RF_SCORE"], lineDict["SpliceAI"]
        ]

        print("\t".join(output))
        results.write("\t".join(output) + "\n")
print("完成annovar结果过滤与格式调整")
