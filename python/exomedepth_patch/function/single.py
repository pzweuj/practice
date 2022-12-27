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
def main(bam, bed, output, reference):
    now = getAbsPath()
    script = now + "/exomeDepthSingle.R"
    cmd = """
        Rscript {script} -b {bed} -i {bam} -o {output} -r {reference}
    """.format(script=script, bed=bed, output=output, reference=reference, bam=bam)
    print(cmd)
    os.system(cmd)

if __name__ == "__main__":
    parser = argparse.ArgumentParser(
        description="""
            Single sample Analysis
            单样本分析，需要导入参考样本rdata
        """,
        prog="single.py",
        usage="python3 single.py [-h] -i <bam file> -b <bed file> -o <output txt> -r <reference rdata>",
        formatter_class=argparse.RawTextHelpFormatter
    )
    parser.add_argument("-v", "--version", action="version",
    version="Version 0.1 20221227")
    parser.add_argument("-i", "--input", type=str,
        help="bam文件")
    parser.add_argument("-b", "--bed", type=str,
        help="bed文件")
    parser.add_argument("-o", "--output", type=str,
        help="结果文件")
    parser.add_argument("-r", "--reference", type=str,
        help="参考文件rdata")
    if len(sys.argv[1:]) == 0:
        parser.print_help()
        parser.exit()
    args = parser.parse_args()
    main(bam=args.input, bed=args.bed, output=args.output, reference=args.reference)

# end
