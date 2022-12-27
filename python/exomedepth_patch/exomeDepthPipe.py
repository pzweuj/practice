#!/usr/bin/env python3
# coding=utf-8
# pzw
# 20221226
# ExomeDepth


import os
import sys
import subprocess
import argparse
import argcomplete

# 获得运行路径
def getAbsPath():
    now = os.path.abspath(os.path.dirname(sys.argv[0]))
    return now

# 主方法
def main(function, option):
    now = getAbsPath()
    function_dict = {
        "single": "single.py",
        "multi": "multiple.py",
        "reference": "reference.py"
    }
    print("python3 " + now + "/function/" + function_dict[function] + " " + " ".join(option))
    subprocess.call("python3 " + now + "/function/" + function_dict[function] + " " + " ".join(option), shell=True)

# 调用
if __name__ == "__main__":
    parser = argparse.ArgumentParser(
        description = "ExomeDepth Pipe",
        prog = "exomeDepthPipe.py",
        usage = "python3 exomeDepthPipe.py [-h] <function> <function option>",
        formatter_class = argparse.RawTextHelpFormatter
    )
    parser.add_argument("-v", "--version", action="version", help="Version 0.1 20221226")
    parser.add_argument("function", choices=("single", "multi", "reference"),
        help="""
            single                单样本模式
            multi                 多样本模式
            reference             基线模式
    """)
    parser.add_argument("option", nargs=argparse.REMAINDER, metavar="function option")
    argcomplete.autocomplete(parser)
    if len(sys.argv[1:]) == 0:
        parser.print_help()
        parser.exit()
    args = parser.parse_args()
    main(args.function, args.option)

# end

