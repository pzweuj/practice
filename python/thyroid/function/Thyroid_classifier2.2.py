# coding=utf-8
# pzw
# 20200526
# 不进行丰度过滤的调整测试
# 增加ROC曲线绘制
# 增加特征输出

import os
import sys
import argparse
from sklearn.ensemble import RandomForestClassifier
from sklearn.metrics import roc_curve, auc
import matplotlib.pyplot as plt
from sklearn.naive_bayes import ComplementNB
import pandas as pd
import joblib
from sklearn.tree import export_graphviz
import numpy as np

## dataframe初始化
def normalization(arff):
	df = pd.read_table(arff, sep="\t", header=0, index_col=False)
	for col in df.columns.values.tolist():
		if col == "Class":
			for i in df["Class"].index:
				if df.loc[i, "Class"] == "pathogenic":
					df.loc[i, "Class"] = 1
				else:
					df.loc[i, "Class"] = 0
	return df

## 重要位点校正。当检测到BRAF V600E，判断为恶性
def BRAFV600E(df, uniqueID):
	BRAF_list = []
	for i in df[uniqueID].index:
		if df.loc[i, uniqueID] == 1:
			BRAF_list.append(df.loc[i, "SampleName"])
	return BRAF_list

## 模型保存
def model_save(clf, model):
	joblib.dump(clf, model)

## 模型读取
def model_load(model):
	clf = joblib.load(model)
	return clf

## 结果预测
def results_predict(clf, x, y, name):
	predict_clf = list(clf.predict(x))
	act_list = list(y)
	sampleName = list(name)
	zipped = list(zip(act_list, predict_clf))
	pred_proba = list(clf.predict_proba(x))
	n_sample = len(zipped)
	ns = 0
	resultsDict = {}
	while ns < n_sample:
		resultsDict[sampleName[ns]] = list(zipped[ns])
		resultsDict[sampleName[ns]].append(list(pred_proba[ns]))
		ns += 1
	return resultsDict

## 模型拟合计算
def model_cal(resultsDict):
	TP = 0
	TN = 0
	FP = 0
	FN = 0
	for i in resultsDict.keys():
		if resultsDict[i][0] == 1:
			if resultsDict[i][1] == 1:
				TP += 1
			else:
				FN += 1
		else:
			if resultsDict[i][1] == 1:
				FP += 1
			else:
				TN += 1
	TPR = "%.3f" % (float(TP) / (TP + FN))
	TNR = "%.3f" % (float(TN) / (TN + FP))
	PPV = "%.3f" % (float(TP) / (TP + FP))
	NPV = "%.3f" % (float(TN) / (TN + FN))
	accuracy = "%.3f" % (float(TP + TN) / (TP + TN + FP + FN))
	print("\t预测恶性\t预测良性")
	print("实际恶性\t" + str(TP) + "\t" + str(FN))
	print("实际良性\t" + str(FP) + "\t" + str(TN) + "\n")
	print("灵敏度TPR\t" + TPR)
	print("特异度TNR\t" + TNR)
	print("PPV\t" + PPV)
	print("NPV\t" + NPV)
	print("准确率\t" + accuracy)

## 位点贡献计算
def feature_ranking(xtrain, clf):
	importances = clf.feature_importances_
	std = np.std([tree.feature_importances_ for tree in clf.estimators_], axis=0)
	indices = np.argsort(importances)[::-1]
	print("Feature ranking:")
	for f in range(xtrain.shape[1]):
		print("%d. feature %d (%f)" % (f + 1, indices[f], importances[indices[f]]))

	plt.figure()
	plt.title("Feature importances")
	plt.bar(range(xtrain.shape[1]), importances[indices],
		color="r", yerr=std[indices], align="center")
	plt.xticks(range(xtrain.shape[1]), indices)
	plt.xlim([-1, xtrain.shape[1]])
	plt.show()


## 随机输出一个决策树
def export_tree(clf, depth, dotFile, pngFile, x):
	export_graphviz(clf.estimators_[depth], out_file=dotFile,
		feature_names=x.columns.values.tolist(),
		class_names=["pathogenic", "benign"], rounded=True,
		proportion=False, precision=2, filled=True
		)
	cmd = """
		dot -Tpng {dotFile} -o {pngFile} -Gdpi=600
	""".format(dotFile=dotFile, pngFile=pngFile)
	os.system(cmd)


## 画ROC
def plot_ROC(test_x, test_y, clf):
	predict_test = clf.predict(test_x)
	prob_test = clf.predict_proba(test_x)
	predict_v = prob_test[:, 1]
	fpr, tpr, _ = roc_curve(test_y, predict_v)
	roc_auc = auc(fpr, tpr)
	plt.title("ROC Curve")
	plt.plot(fpr, tpr, "b", label="AUC = %0.2f" % roc_auc)
	plt.legend(loc="lower right")
	plt.plot([0, 1], [0, 1], "r--")
	plt.xlim([0, 1])
	plt.ylim([0, 1])
	plt.ylabel("True Positive Rate")
	plt.xlabel("False Positive Rate")
	plt.show()

## 主方法
def main(arff, modelFile, resultsFile, method, n_est, depth, classifier, V600E, Tree, ROC, feature):
	
	# 是否训练
	if method == "train":
		## 模型建立，模型设置树数100，深度9，可继续调整获得最佳树数和深度
		if classifier == "Complement_Naive_Bayes":
			model = ComplementNB()
		else:
			model = RandomForestClassifier(n_estimators=n_est, bootstrap=True, max_depth=depth)
		## 读入训练组
		train = normalization(arff)
		xtrain = train.iloc[:, 2:]
		ytrain = train["Class"]
		ytrain = ytrain.astype("int")
		name_train = train["SampleName"]
		## 模型训练
		clf = model.fit(xtrain, ytrain)
		## 模型拟合计算
		resultsDict = results_predict(clf, xtrain, ytrain, name_train)
		model_cal(resultsDict)
		## 模型保存
		model_save(clf, modelFile)

		if feature:
			feature_ranking(xtrain, clf)

	## 测试与预测
	else:
		## 读取模型
		clf = model_load(modelFile)
		## 读入测试组
		test = normalization(arff)
		xtest = test.iloc[:, 2:]
		ytest = test["Class"]
		name_test = test["SampleName"]
		## 结果统计
		resultsDict = results_predict(clf, xtest, ytest, name_test)
		results = open(resultsFile, "w")
		results.write("样本编号\t临床结果\t预测结果\t良性概率\t恶性概率\n")
		
		if V600E != "":
			V600E_list = BRAFV600E(test, V600E)

		for i in resultsDict.keys():
			if resultsDict[i][0] == 1:
				act_result = "pathogenic"
			else:
				act_result = "benign"
			if resultsDict[i][1] == 1:
				predict_result = "pathogenic"
			else:
				predict_result = "benign"

			pred_0 = "%.2f" % resultsDict[i][2][0]
			pred_1 = "%.2f" % resultsDict[i][2][1]

			if V600E != "":
				if i in V600E_list:
					if resultsDict[i][1] == 0:
						resultsDict[i][1] = 1
						predict_result = "pathogenic"
						pred_0 = "0"
						pred_1 = "1"

			output = "\t".join([i, act_result, predict_result, pred_0, pred_1]) + "\n"
			results.write(output)
		results.close()
		model_cal(resultsDict)
		if Tree:
			export_tree(clf, depth, "test.dot", "test.png", xtest)

		if ROC:
			plot_ROC(xtest, ytest, clf)

if __name__ == "__main__":
	parser = argparse.ArgumentParser(
		description="Thyroid Nodule Classification @PZW sklearn版本",
		prog="Thyroid_classifier.py",
		usage="""python3 Thyroid_classifier.py [-h] -a <arffFile> -m <model> [-o <resultsFile>] [-md train]
			\n训练模式使用：
				python3 Thyroid_classifier.py -a <arff> -m <model> -md train
			\n测试模式使用：
				python3 Thyroid_classifier.py -a <arff> -m <model> -o <resultsFile>

		""",
		formatter_class=argparse.RawTextHelpFormatter
	)
	parser.add_argument("-v", "--version", action="version",
		version="Version 0.1 20200519")
	parser.add_argument("-a", "--arff", type=str,
		help="输入训练/测试文件")
	parser.add_argument("-m", "--model", type=str,
		help="输出/输入模型")
	parser.add_argument("-o", "--output", type=str,
		help="预测模式下可选，输出预测结果", default="")
	parser.add_argument("-md", "--method", type=str,
		help="设置'-md train' 将使用训练模式", default="no")
	parser.add_argument("-n", "--n_est", type=int,
		help="可选，设置树数，仅在训练模式下有效，默认为100", default=100)
	parser.add_argument("-d", "--depth", type=int,
		help="可选，设置决策树深度，仅在训练模式下有效，默认为12", default=12)
	parser.add_argument("-c", "--classifier", type=str,
		help="可选，选择训练算法，默认为随机森林。可选'Complement_Naive_Bayes'", default="")
	parser.add_argument("-fix", "--fix", type=str,
		help="BRAF V600E校正", default="")
	parser.add_argument("-tree", "--tree", type=bool,
		help="随机输出一个决策树，默认不输出，只有预测模式可用", default=False)
	parser.add_argument("-ROC", "--ROC", type=bool,
		help="输出ROC曲线，默认不输出，只有预测模式可用", default=False)
	parser.add_argument("-feature", "--feature", type=bool,
		help="输出训练特征，默认不输出，只有训练模式可用", default=False)
	if len(sys.argv[1:]) == 0:
		parser.print_help()
		parser.exit()
	args = parser.parse_args()
	main(arff=args.arff, modelFile=args.model, resultsFile=args.output, method=args.method, n_est=args.n_est, depth=args.depth, classifier=args.classifier, V600E=args.fix, Tree=args.tree, ROC=args.ROC, feature=args.feature)
