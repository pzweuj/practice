# coding=utf-8
# pzw
# 20200304
# Thyroid Nodule Classification
# 用于进行甲状腺结节良恶性判断
# 基于python3

"""
sudo apt-get install python3-pil python3-matplotlib python3-pygraphviz
sudo pip3 install javabridge
sudo pip3 install python-weka-wrapper3
"""

import os
import sys
import argparse
import weka.core.jvm as jvm
from weka.core.converters import Loader
from weka.classifiers import Classifier
from weka.classifiers import FilteredClassifier
from weka.classifiers import Evaluation
from weka.core.classes import Random
import weka.core.serialization as serialization

# vcf变异信息解析
def vcf_data():
	pass

# 超声得分导入
def echo_data():
	pass

def getAbsPath():
	now = os.path.abspath(os.path.dirname(sys.argv[0]))
	return now


# Training Model 使用训练集生成训练矩阵
def TrainingModel(arff, modelOutput, clsfier):
	# 启动java虚拟机
	jvm.start()
	# 导入训练集
	loader = Loader(classname="weka.core.converters.ArffLoader")
	train = loader.load_file(arff)
	train.class_is_first()
	# 使用RandomForest算法进行训练，因为在GUI版本weka中使用多种方式训练后发现此方式TPR与TNR较高
	cls_name = "weka.classifiers." + clsfier
	clsf = Classifier(classname=cls_name)
	clsf.build_classifier(train)
	print(clsf)
	# 建立模型
	fc = FilteredClassifier()
	fc.classifier = clsf
	evl = Evaluation(train)
	evl.crossvalidate_model(fc, train, 10, Random(1))
	print(evl.percent_correct)
	print(evl.summary())
	print(evl.class_details())
	print(evl.matrix())
	# 结果统计
	matrixResults = evl.confusion_matrix
	TN = float(matrixResults[0][0])
	FP = float(matrixResults[0][1])
	FN = float(matrixResults[1][0])
	TP = float(matrixResults[1][1])
	TPR = TP / (TP + FN)
	TNR = TN / (FP + TN)
	PPV = TP / (TP + FP)
	NPV = TN / (TN + FN)
	print("算法： " + clsfier)
	print("敏感度 TPR: " + str(TPR))
	print("特异度 TNR: " + str(TNR))
	print("PPV: " + str(PPV))
	print("NPV: " + str(NPV))
	# 保存模型
	clsf.serialize(modelOutput, header=train)
	# 退出虚拟机
	jvm.stop()
	print("分析模型建立完成")

def TestClassification(arff, modelInput, results, sampleName):
	# 启动java虚拟机
	jvm.start()
	# 导入分析模型
	objects = serialization.read_all(modelInput)
	clsf = Classifier(jobject=objects[0])
	print(clsf)
	# 导入测试组
	loader = Loader(classname="weka.core.converters.ArffLoader")
	test = loader.load_file(arff)
	test.class_is_first()
	# 分析结果
	resultsFile = open(results, "w")
	if sampleName:
		resultsFile.write("样本编号\t原判断\t预测\t良性概率\t恶性概率\n")
		print("样本编号\t原判断\t预测\t良性概率\t恶性概率")
		sampleNameListFile = open(getAbsPath() + "/temp.txt", "r")
		sampleNameList = []
		for snlf in sampleNameListFile:
			sampleNameList.append(snlf.split("\n")[0])
		sampleNameListFile.close()
	else:
		resultsFile.write("序号\t原判断\t预测\t良性概率\t恶性概率\n")
		print("序号\t原判断\t预测\t良性概率\t恶性概率")
	for index, inst in enumerate(test):
		pred = clsf.classify_instance(inst)
		dist = clsf.distribution_for_instance(inst)
		if sampleName:
			sampleID = sampleNameList[index]
		else:
			sampleID = str(index + 1)
		origin = inst.get_string_value(inst.class_index)
		prediction = inst.class_attribute.value(int(pred))
		sameAsOrigin = "yes" if pred != inst.get_value(inst.class_index) else "no"
		NRate = dist.tolist()[0]
		PRate = dist.tolist()[1]
		resultsFile.write("%s\t%s\t%s\t%s\t%s" % (sampleID, origin, prediction, str(NRate), str(PRate)) + "\n")
		print("%s\t%s\t%s\t%s\t%s" % (sampleID, origin, prediction, str(NRate), str(PRate)))
	resultsFile.close()
	# 退出java虚拟机
	jvm.stop()
	print("检测完成")

def main(arff, modelFile, results, method, clsfier, sampleName):
	if method == "train":
		TrainingModel(arff, modelFile, clsfier)
	else:
		print("if use train sets, please specify '-md train'")
		TestClassification(arff, modelFile, results, sampleName)

if __name__ == "__main__":
	parser = argparse.ArgumentParser(
		description="Thyroid Nodule Classification @PZW",
		prog="thyroid_nodule_classification.py",
		usage="python3 thyroid_nodule_classification.py [-h] -a <arffFile> -m <model> [-o <resultsFile>] [-md train]",
		formatter_class=argparse.RawTextHelpFormatter
	)
	parser.add_argument("-v", "--version", action="version",
		version="Version 0.1 20200304")
	parser.add_argument("-a", "--arff", type=str,
		help="input the arff format file")
	parser.add_argument("-m", "--model", type=str,
		help="input/output a model file")
	parser.add_argument("-o", "--output", type=str,
		help="test sets prediction output", default="")
	parser.add_argument("-md", "--method", type=str,
		help="'-md train' switch to train process", default="no")
	parser.add_argument("-c", "--classifier", type=str,
		help="specify a classifier, default=trees.RandomForest", default="trees.RandomForest")
	parser.add_argument("-s", "--sampleName", type=bool,
		help="output sample name, '-s True'", default=False)

	if len(sys.argv[1:]) == 0:
		parser.print_help()
		parser.exit()
	args = parser.parse_args()
	main(arff=args.arff, modelFile=args.model, results=args.output, method=args.method, clsfier=args.classifier, sampleName=args.sampleName)