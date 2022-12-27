# coding=utf-8
# pzw
# 20221226

import os
import sys
import argparse

# 获得运行路径
def getAbsPath():
    now = os.path.abspath(os.path.dirname(sys.argv[0]))
    return now

# 主方法
def main(bamDir, bed, outputDir):
    now = getAbsPath()
    script = now + "/exomeDepthPipe.R"
    cmd = """
        Rscript {script} -i {bamDir} -b {bed} -o {outputDir}
    """.format(script=script, bamDir=bamDir, bed=bed, outputDir=outputDir)
    print(cmd)
    os.system(cmd)

if __name__ == "__main__":
    parser = argparse.ArgumentParser(
        description="""
            Multiple samples Analysis
            多样本分析，每个样本互为参考
        """,
        prog="multiple.py",
        usage="python3 multiple.py [-h] -i <bams directory> -b <bed file> -o <output directory>",
        formatter_class=argparse.RawTextHelpFormatter
    )
    parser.add_argument("-v", "--version", action="version",
    version="Version 0.1 20221226")
    parser.add_argument("-i", "--input", type=str,
        help="bam文件存放文件夹")
    parser.add_argument("-b", "--bed", type=str,
        help="bed文件")
    parser.add_argument("-o", "--output", type=str,
        help="结果输出文件夹")
    if len(sys.argv[1:]) == 0:
        parser.print_help()
        parser.exit()
    args = parser.parse_args()
    main(bamDir=args.input, bed=args.bed, outputDir=args.output)

# end
