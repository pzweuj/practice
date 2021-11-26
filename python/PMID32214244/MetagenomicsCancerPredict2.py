# coding=utf-8
# pzw
# 20211125

import pandas as pd
from math import log
import joblib
from sklearn.preprocessing import StandardScaler
from sklearn.svm import LinearSVC
from sklearn.pipeline import make_pipeline
import seaborn as sns
import matplotlib.pyplot as plt

## 模型保存
def model_save(clf, model):
	joblib.dump(clf, model)

## 模型读取
def model_load(model):
	clf = joblib.load(model)
	return clf

# 分类癌症与正常
def normalization(table, clsList):
    df = pd.read_table(table, sep="\t", header=0, index_col=0)
    for i in df["Class"].index:
        types = df.loc[i, "Class"]
        for c in range(len(clsList)):
            if types in clsList[c].split(","):
                df.loc[i, "Class"] = c
    
    dropRow = []
    for i in df.index:
        if not isinstance(df.loc[i, "Class"], int):
            dropRow.append(i)
    df = df.drop(dropRow)

    colList = df.columns.values.tolist()
    colSumDict = {}
    for i in df.index:
        gSum = 0
        for col in colList:
            if col != "Class" and col != "分组":
                gSum += df.loc[i, col]
        colSumDict[i] = gSum

    for i in df.index:
        for col in colList:
            if col != "Class" and col != "分组":
                # logCPM = log2(CPM + 1)
                # https://support.bioconductor.org/p/107719/
                # 可能会导致倍数错误
                df.loc[i, col] = log((df.loc[i, col] / colSumDict[i] * 1000000 + 1), 2)

    return df

# 模型测试
def testModel(xtest, clf, clsList):
    predict_clf = list(clf.predict(xtest))
    predict_transform = []
    for pc in predict_clf:
        predict_transform.append(clsList[pc])
    name_test = xtest.index
    sampleName = list(name_test)
    zipped = list(zip(sampleName, predict_transform))
    return zipped

# 贡献度统计
def feature_importance(x, model, clsName):
    feature_names = x.columns.values.tolist()
    coefs = model.named_steps[clsName].coef_.flatten()
    zipped = zip(feature_names, coefs)
    df = pd.DataFrame(zipped, columns=["feature", "value"])
    df["abs_value"] = df["value"].apply(lambda x: abs(x))
    df["colors"] = df["value"].apply(lambda x: "green" if x > 0 else "red")
    df = df.sort_values("abs_value", ascending=False)
    return df

# 贡献度作图
def feature_plot(df, picName):
    fig, ax = plt.subplots(1, 1, figsize=(15, 12))
    sns.barplot(x="feature", y="value", data=df.head(20), palette=df.head(20)["colors"])
    ax.set_xticklabels(ax.get_xticklabels(), rotation=90, fontsize=20)
    ax.set_title("Top 20 Features", fontsize=25)
    ax.set_ylabel("Coef", fontsize=22)
    ax.set_xlabel("Feature Name", fontsize=20)
    plt.savefig(picName)

# 模型训练
def ownLinearSVC(xtrain, ytrain):
    model = make_pipeline(StandardScaler(), LinearSVC(random_state=0, tol=1e-4, max_iter=1000, C=1.0))
    clf = model.fit(xtrain, ytrain)
    return clf

# 数据分析
####################################################
trainTable = "CountMatrix/TCGA.train.txt"
testTable = "CountMatrix/TCGA.test.txt"
saveModel = "Model/TCGA.LinearSVC.OtherVsSKCM.model"
results = "Results/TCGA.LinearSVC.OtherVsSKCM.results.txt"
featureFile = "FeatureRank/TCGA.LinearSVC.OtherVsSKCM.FeatureRank.txt"
featurePic = "FeatureRank/TCGA.LinearSVC.OtherVsSKCM.FeatureRank.Top20.png"
typeList = ["Control,PRAD,NSCLC", "SKCM"]
####################################################

## 分析流程
train = normalization(trainTable, typeList)
print(train)
xtrain = train.iloc[:, 2:]
ytrain = train["Class"]
ytrain = ytrain.astype("int")
test = normalization(testTable, typeList)
xtest = test.iloc[:, 2:]
ytest = test["Class"]
ytest = ytest.astype("int")
clf = ownLinearSVC(xtrain, ytrain)
df = feature_importance(xtrain, clf, "linearsvc")
print(df.head(20))
df.to_csv(featureFile, sep="\t", index=False)
feature_plot(df, featurePic)
model_save(clf, saveModel)
clf = model_load(saveModel)
zList = testModel(xtest, saveModel, typeList)
r = open(results, "w", encoding="utf-8")
for z in zList:
    r.write(z[0] + "\t" + str(z[1]) + "\n")
r.close()

