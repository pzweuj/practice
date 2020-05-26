# coding=utf-8
# pzw
# 20200330

import sys
import subprocess
import argparse
import argcomplete
import os


def main(function, option):

	now = os.path.abspath(os.path.dirname(sys.argv[0]))
	function_dict = {
		"call_variants": "/function/thyroid_variants_call.py",
		"make_variants_matrix": "/function/thyroid_matrix.py",
		"matrix_filter": "/function/thyroid_matrix_filter.py",
		"arff_maker": "/function/thyroid_data3.py",
		"classification": "/function/Thyroid_classifier2.1.py"
	}
	print("python3 " + now + function_dict[function] + " " + " ".join(option))
	subprocess.call("python3 " + now + function_dict[function] + " " + " ".join(option), shell=True)

if __name__ == "__main__":
	now=os.path.dirname(sys.executable)
	parser = argparse.ArgumentParser(description="Thyroid Function", prog="ThyroidPro.py",
		usage="python3 ThyroidPro.py [-h] <function> <function option>", formatter_class=argparse.RawTextHelpFormatter)
	parser.add_argument("-v", "--version", action="version", version="Version 0.1")
	parser.add_argument("function",
		choices=(
			"call_variants",
			"make_variants_matrix",
			"matrix_filter",
			"arff_maker",
			"classification"
			),
		help="""
\033[32;1m----------------甲状腺结节分析模块----------------\033[0m
call_variants                   变异检测
make_variants_matrix            生成变异位点矩阵
matrix_filter                   变异位点矩阵过滤（建议过滤后再进行一次手动过滤）
arff_maker                      生成arff文件
classification                  机器学习分类器

""")
	parser.add_argument("option", nargs=argparse.REMAINDER, metavar="function option")
	argcomplete.autocomplete(parser)
	if len(sys.argv[1:]) == 0:
		parser.print_help()
		parser.exit()
	args = parser.parse_args()
	main(function=args.function, option=args.option)
