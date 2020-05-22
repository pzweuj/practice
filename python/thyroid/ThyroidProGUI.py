# coding=utf-8
# pzw
# 20200520
# 甲状腺结节良恶性分析
# GUI版本
# v0.2

import PySimpleGUI as sg
import os
# import ThyroidPro
import function.thyroid_variants_call
import function.thyroid_matrix
import function.thyroid_matrix_filter
import function.thyroid_data3
import function.Thyroid_classifier2

# 改变窗体主题颜色
sg.ChangeLookAndFeel("GreenTan")

# 变异检测参数流程
muts_detect_layout = [
	[sg.Text("甲状腺数据变异检测")],
	[sg.Text("输入样本列表"), sg.Input(key="sampleList"), sg.FileBrowse()],
	[sg.Text("原始数据文件夹"), sg.Input(key="rawdataDir"), sg.FolderBrowse()],
	[sg.Text("输出文件夹"), sg.Input(key="outputDir"), sg.FolderBrowse()],
	[sg.Text("线程数"), sg.Combo(["1", "2", "3", "4", "5", "6", "7", "8"], key="threads", default_value="8")],
	[sg.Submit()]
	# [sg.Text("Log"), sg.Output(key="muts_log")]
]

# 矩阵参数流程
matrix_layout = [
	[sg.Text("变异矩阵生成")],
	[sg.Text("输入注释数据文件夹"), sg.Input(key="annovarDir"), sg.FolderBrowse()],
	[sg.Text("输出矩阵文件"), sg.Input(key="matrixFile"), sg.FileBrowse()],
	[sg.Submit()]
	# [sg.Text("Log"), sg.Output(key="matrix_log")]
]


# 变异矩阵过滤
clean_matrix_layout = [
	[sg.Text("变异矩阵过滤")],
	[sg.Text("输入矩阵文件"), sg.Input(key="matrixFile2"), sg.FileBrowse()],
	[sg.Text("输出过滤后矩阵"), sg.Input(key="cleanMatrixFile"), sg.FileBrowse()],
	[sg.Submit()]
	# [sg.Text("Log"), sg.Output(key="clean_log")]
]

# Arff文件生成
arff_layout = [
	[sg.Text("arff文件生成")],
	[sg.Text("注释结果文件夹"), sg.Input(key="annovarDir2"), sg.FolderBrowse()],
	[sg.Text("输出arff文件"), sg.Input(key="arff"), sg.FileBrowse()],
	[sg.Text("输入bed文件"), sg.Input(key="bed"), sg.FileBrowse()],
	[sg.Text("良恶性信息文件"), sg.Input(key="info"), sg.FileBrowse()],
	[sg.Text("注：测试arff文件生成可以不加入良恶性信息！")],
	[sg.Submit()]
]

# 机器学习训练器
train_layout = [
	[sg.Text("机器学习训练器")],
	[sg.Text("输入arff文件"), sg.Input(key="arff2"), sg.FileBrowse()],
	[sg.Text("输出分析模型"), sg.Input(key="model"), sg.FileBrowse()],
	[
		sg.Text("选择训练算法"),
		sg.Combo([
				"RandomForest",
				"Complement_Naive_Bayes"
			], key="clsMethod", default_value="RandomForest")
	],
	[sg.Submit()]
	# [sg.Text("Log"), sg.Output(key="train_log")]
]


# 机器学习分类器
test_layout = [
	[sg.Text("机器学习预测器")],
	[sg.Text("输入arff文件"), sg.Input(key="arff3"), sg.FileBrowse()],
	[sg.Text("输入训练模型"), sg.Input(key="model2"), sg.FileBrowse()],
	[sg.Text("输出预测结果"), sg.Input(key="results"), sg.FileBrowse()],
	[sg.Text("BRAF V600E校正点"), sg.Input(key="V600E")],
	[sg.Submit()]
	# [sg.Text("Log"), sg.Output(key="test_log")]
]

# 总框架
layout = [
	[sg.Text("甲状腺结节良恶性判断分析")],
	[sg.TabGroup([
		[sg.Tab("变异检测", muts_detect_layout)],
		[sg.Tab("矩阵生成", matrix_layout)],
		[sg.Tab("矩阵过滤", clean_matrix_layout)],
		[sg.Tab("Arff生成", arff_layout)],
		[sg.Tab("训练", train_layout)],
		[sg.Tab("预测", test_layout)]
		])
	],
	[sg.Cancel("Exit")]
]


window = sg.Window("甲状腺结节良恶性判断分析", layout)
stuff = window.read()

event = stuff[0]
runningDict = stuff[1]

if event in (None, "Exit"):
	window.Close()
else:
	tag = runningDict[0]
	if tag == "变异检测":
		function.thyroid_variants_call.main(runningDict["rawdataDir"], runningDict["outputDir"], runningDict["sampleList"], runningDict["threads"])
	elif tag == "矩阵生成":
		function.thyroid_matrix.main(runningDict["annovarDir"], runningDict["matrixFile"])
	elif tag == "矩阵过滤":
		function.thyroid_matrix_filter.main(runningDict["matrixFile2"], runningDict["cleanMatrixFile"])
	elif tag == "Arff生成":
		function.thyroid_data3.main(runningDict["annovarDir2"], runningDict["arff"], runningDict["bed"], runningDict["info"])
	elif tag == "训练":
		function.Thyroid_classifier2.main(runningDict["arff2"], runningDict["model"], "", "train", 100, 12, runningDict["clsMethod"], "")
	elif tag == "预测":
		function.Thyroid_classifier2.main(runningDict["arff3"], runningDict["model2"], runningDict["results"], "", 100, 12, runningDict["clsMethod"], runningDict["V600E"])
	else:
		windows.Close()
