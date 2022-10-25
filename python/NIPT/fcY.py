# coding=utf-8
# pzw
# 20221021
# Y染色体浓度计算
# 使用loess校正后数据计算Y浓度
# Chiu et al, Non-invasive prenatal assessment of trisomy 21 by multiplexed maternal plasma DNA sequencing: large scale validity study, 2011

import pandas as pd
import os

# 读入结果并提取需要的列
def GetURDateFrame(corFile):
    df = pd.read_csv(corFile, header=None, sep="\t")
    # df.rename(columns={0: "chrom", 1: "start", 2: "end", 3: "name", 5: "GC", 13: "UR", 15: "URGC"}, inplace=True)
    df.rename(columns={0: "chrom", 1: "start", 2: "end", 4: "GC", 12: "UR", 14: "URGC"}, inplace=True)
    # df2 = df[["chrom", "start", "end", "name", "GC", "UR", "URGC"]]
    df2 = df[["chrom", "start", "end", "GC", "UR", "URGC"]]

    # UR均值
    URmean = df2["UR"].mean()

    # UR修正值
    df2.insert(loc=6, column="URAdj", value=0)
    nonZeroList = df2.loc[df2["URGC"] != 0].index
    for i in nonZeroList:
        URGC = df2.loc[i, "URGC"]
        UR = df2.loc[i, "UR"]
        URAdj = round(UR * (URmean / URGC))
        df2.loc[i, "URAdj"] = URAdj
    return df2

# 获得Y丰度
def percY(dataFrame):
    sumY = dataFrame[dataFrame["chrom"] == "chrY"]["URAdj"].sum()
    sumAll = dataFrame["URAdj"].sum()
    # sumParX = dataFrame[dataFrame["name"].str.contains("chromX.PAR")]["URAdj"].sum()
    # sumParY = dataFrame[dataFrame["name"].str.contains("chromY.PAR")]["URAdj"].sum()
    Y_per = float(sumY) / sumAll
    return Y_per

# 获得一批次的平均Y丰度
def meanPerY(dir):
    Y_per_List = []
    for i in os.listdir(dir):
        df = GetURDateFrame(dir + "/" + i)
        y_p = percY(df)
        Y_per_List.append(y_p)
    meanYPer = sum(Y_per_List) / len(Y_per_List)
    return meanYPer



# Main
Y_per_female_mean = 0.0016487657371067291
Y_per_male_mean = 0.0022046226182481808


# 批量
output = open("ffy_output.txt", "w", encoding="utf-8")
for i in os.listdir("corCounts/train_male"):
    name = i.split(".")[0]
    df = GetURDateFrame("corCounts/train_male/" + i)
    s = percY(df)
    ffy = (s - Y_per_female_mean) / (Y_per_male_mean - Y_per_female_mean)
    


