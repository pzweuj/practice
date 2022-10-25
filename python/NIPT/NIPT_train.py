# coding=utf-8
# pzw
# 20221021
# seqFF
# https://pubmed.ncbi.nlm.nih.gov/25967380/
# 预测胎儿性别
# 使用弹性网络+分级加权

import os
import pandas as pd
from sklearn.linear_model import ElasticNetCV
import joblib

# 模型保存
def model_save(clf, model):
    joblib.dump(clf, model)

# 模型读取
def model_load(model):
    clf = joblib.load(model)
    return clf


# 读取多个样本，形成一个大matrix
def corExact(corFile, ffDict):
    fi = corFile
    name = fi.split("/")[-1].split(".")[0]
    ff = ffDict[name]
    df = pd.read_csv(fi, header=None, sep="\t")
    df.rename(columns={0: "chrom", 1: "start", 2: "end", 4: "GC", 12: "UR", 14: "URGC"}, inplace=True)
    df2 = df[["chrom", "start", "end", "GC", "UR", "URGC"]]
    URmean = df2["UR"].mean()
    df2.insert(loc=6, column="URAdj", value=0)
    nonZeroList = df2.loc[df2["URGC"] != 0].index
    for i in nonZeroList:
        URGC = df2.loc[i, "URGC"]
        UR = df2.loc[i, "UR"]
        URAdj = round(UR * (URmean / URGC))
        df2.loc[i, "URAdj"] = URAdj
    df3 = df2[(df2["chrom"] != "chrY") & (df2["chrom"] != "chrX") & (df2["chrom"] != "chr13") & (df2["chrom"] != "chr18") & (df2["chrom"] != "chr21")]
    df3["checkPoint"] = df3["chrom"] + "-" + df3["start"].astype(str)
    df4 = df3[["checkPoint", "URAdj"]]
    df5 = df4.reset_index(drop=True)
    lastIndex = len(df5.index)
    df5.loc[lastIndex, "checkPoint"] = "ff"
    df5.loc[lastIndex, "URAdj"] = ff
    df5.rename(columns={"URAdj": name}, inplace=True)
    return df5


# 创建矩阵
def InputMatrix(dir, ffDict):
    n = 0
    for i in os.listdir(dir):
        n += 1
        f = dir + "/" + i
        if n == 1:
            d1 = corExact(f, ffDict)
        else:
            di = corExact(f, ffDict)
            mergeDf = pd.merge(d1, di, on="checkPoint")
            d1 = mergeDf
    return d1



############## 输出counts矩阵 ############
# 读取胎儿浓度
# ffDict = {}
# with open("all_sample_ff.txt") as f:
#     for line in f:
#         if not line.startswith("#"):
#             lines = line.split("\t")
#             ffDict[lines[0]] = float(lines[1])

# df = InputMatrix("corCounts/train_male", ffDict)
# df.to_csv("output.txt", sep="\t", header=True, index=False)
########################################



################ 训练 ##################
# seqFF训练用的样本数：25312
# seqFF中的两个常数
# parameter1: -0.000315420152542636
# parameter2: 129.062773589499

################# 训练 ####################
df = pd.read_csv("train_male.txt", sep="\t", header=0, index_col=0)
dft = df.T
cols = [i for i in dft.columns if i != "ff"]
X = dft[cols]
y = dft["ff"]
regr = ElasticNetCV(cv=10, random_state=0, max_iter=1000)
regr.fit(X, y)
model_save(regr, "train_male.model")

################# 测试 #####################
df_test = pd.read_csv("test_female.txt", sep="\t", header=0, index_col=0)
dft_test = df_test.T
cols_test = [i for i in dft_test.columns if i != "ff"]
X_test = dft_test[cols]
y_test_act = list(dft_test["ff"])
regr = model_load("train_male.model")
predict_test = list(regr.predict(X_test))
zipped = list(zip(y_test_act, predict_test))
for z in zipped:
    print(z[0], z[1])


