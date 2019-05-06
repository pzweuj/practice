# coding=utf-8
# pzw
# 20190506
# v0.1

import sys
import os
import argparse
import subprocess
import argcomplete

def main(function, option):

	now = os.path.dirname(sys.executable)

	function_dict = {
		"reference": "function/BRCA_reference_v1.0.py",
		"calculation": "function/BRCA_analysis_v1.0.py",
		"summary": "function/BRCA_Long_Del_v0.1.py",
		"plot": "function/BRCA_plot_v0.1.py"
	}

	print("python3 " + os.path.join(now, function_dict[function]) + " " + " ".join(option))
	subprocess.call("python3 " + os.path.join(now, function_dict[function]) + " " + " ".join(option), shell=True)

if __name__ == "__main__":
	now = os.path.dirname(sys.executable)
	parser = argparse.ArgumentParser(description="BRCA Long Deletions Analysis Toolkit",
		prog="BRCA_toolkit_v0.1.py",
		usage="python3 BRCA_toolkit_v0.1.py [-h] <function> <function option>",
		formatter_class=argparse.RawTextHelpFormatter
	)
	parser.add_argument("-v", "--version", action="version", version="Version 0.1 20190506")
	parser.add_argument("function", choices=("reference", "calculation", "summary", "plot"),
		help='''
			reference                    创建参考文件
			calculation                  计算Z值或DQ值
			summary                      得到最终结果
			plot                         可视化
	''')
	parser.add_argument("option", nargs=argparse.REMAINDER, metavar="function option")
	argcomplete.autocomplete(parser)
	args = parser.parse_args()
	main(function=args.function, option=args.option)


