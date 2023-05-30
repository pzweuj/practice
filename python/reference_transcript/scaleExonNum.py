# coding=utf-8
# pzw
# 20230523
# 从UCSC下载的外显子bed中，修正外显子的编号，补充内含子
# 补充基因名，仅保留参考转录本

# import os
import sys

# 结构是 {"转录本": "基因"}
def ref_trans(refFile):
    refDict = {}
    with open(refFile, "r", encoding="utf-8") as f:
        for line in f:
            if not line.startswith("#"):
                lines = line.split("\t")
                refDict[lines[1].split(".")[0]] = lines[0]
    return refDict

# 处理转录本文件
def transcript_store(ucsc_exon_bed, refFile):
    transcript_dict = ref_trans(refFile)
    # print("【Step1】参考转录本读取完成！")
    ref_trans_list = list(transcript_dict.keys())
    rawBed = open(ucsc_exon_bed, "r", encoding="utf-8")
    location_dict = {}
    for line in rawBed:
        transcript_check = line.split("\t")[3].split("_exon")[0].split(".")[0]
        if not transcript_check in ref_trans_list:
            # print(transcript_check, "非参考转录本，跳过")
            continue
        else:
            # print(transcript_check, "是参考转录本，录入信息中")
            lines = line.replace("\n", "").split("\t")
            chrom = lines[0]

            # 过滤掉多态性的小染色体
            if "_" in chrom:
                continue

            start = lines[1]
            end = lines[2]
            transcript = lines[3].split("_exon")[0]
            if not transcript in location_dict.keys():
                gene = transcript_dict[transcript.split(".")[0]]
                strand = lines[5]
                location_dict[transcript] = {"start": [start], "end": [end], "strand": strand, "chrom": chrom, "gene": gene}
            else:
                location_dict[transcript]["start"].append(start)
                location_dict[transcript]["end"].append(end)
    rawBed.close()
    # print("【Step2】参考转录本坐标录入完成！")
    return location_dict

# 补充内含子
def intron_append(location_dict):
    for k in location_dict:
        if not "intron_start" in location_dict[k]:
            location_dict[k]["intron_start"] = []
        if not "intron_end" in location_dict[k]:
            location_dict[k]["intron_end"] = []
        zipped = list(zip(location_dict[k]["start"], location_dict[k]["end"]))
        exonNum = len(zipped)
        
        for i in zipped:
            # 每个外显子的终止坐标+1就是内含子的起始坐标
            location_dict[k]["intron_start"].append(str(int(i[1]) + 1))
            # 每个外显子的起始坐标-1就是内含子的终止坐标
            location_dict[k]["intron_end"].append(str(int(i[0]) - 1))
        
        # 然后把内含子最后一个起始坐标和第一个终止坐标舍弃
        location_dict[k]["intron_start"] = location_dict[k]["intron_start"][:-1]
        location_dict[k]["intron_end"] = location_dict[k]["intron_end"][1:]
    # print("【Step3】内含子坐标补充完成！")
    return location_dict

# 最后根据转录方向修改坐标的方向
def strand_repair(location_dict):
    for k in location_dict:
        strand = location_dict[k]["strand"]
        if strand == "-":
            location_dict[k]["start"] = location_dict[k]["start"][::-1]
            location_dict[k]["end"] = location_dict[k]["end"][::-1]
            location_dict[k]["intron_start"] = location_dict[k]["intron_start"][::-1]
            location_dict[k]["intron_end"] = location_dict[k]["intron_end"][::-1]
    # print("【Step4】转录方向校对完成！")
    return location_dict

# 写入文件
def write_bed(location_dict, outputName):
    output = open(outputName, "w", encoding="utf-8")
    for k in location_dict:
        gene = location_dict[k]["gene"]
        strand = location_dict[k]["strand"]
        start = location_dict[k]["start"]
        end = location_dict[k]["end"]
        chrom = location_dict[k]["chrom"]
        intron_start = location_dict[k]["intron_start"]
        intron_end = location_dict[k]["intron_end"]
        for s in range(len(start)):
            outputStringList = [chrom, start[s], end[s], "exon" + str(s+1), gene, k, strand]
            output.write("\t".join(outputStringList) + "\n")
        for i in range(len(intron_start)):
            outputStringList = [chrom, intron_start[i], intron_end[i], "intron" + str(i+1), gene, k, strand]
            output.write("\t".join(outputStringList) + "\n")
    output.close()
    # print("【Step5】写入文件完成！")

# main
def main(inputUcscFile, outputFileName, refTransFile):
    transcript_count_dict = transcript_store(inputUcscFile, refTransFile)
    intron_dict = intron_append(transcript_count_dict)
    strand_dict = strand_repair(intron_dict)
    write_bed(strand_dict, outputFileName)

########################################
try:
    main(sys.argv[1], sys.argv[2], sys.argv[3])
except:
    print("python scaleExonNum.py <ucsc exon bed> <output exon/intron bed> <ref transcript file>")

# end
