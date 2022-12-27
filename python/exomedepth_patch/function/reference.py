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
def main(bamDir, bed, output):
    now = getAbsPath()
    script = now + "/exomeDepthRef.R"
    cmd = """
        Rscript {script} -i {bamDir} -b {bed} -o {output}
    """.format(script=script, bamDir=bamDir, bed=bed, output=output)
    os.system(cmd)

if __name__ == "__main__":
    parser = argparse.ArgumentParser(
        description="""
            Reference create
            参考样本count保存
        """,
        prog="reference.py",
        usage="python3 reference.py [-h] -i <bams directory> -b <bed file> -o <output rdata>",
        formatter_class=argparse.RawTextHelpFormatter
    )
    parser.add_argument("-v", "--version", action="version",
    version="Version 0.1 20221226")
    parser.add_argument("-i", "--input", type=str,
        help="bam文件存放文件夹")
    parser.add_argument("-b", "--bed", type=str,
        help="bed文件")
    parser.add_argument("-o", "--output", type=str,
        help="结果rdata生成")
    if len(sys.argv[1:]) == 0:
        parser.print_help()
        parser.exit()
    args = parser.parse_args()
    main(bamDir=args.input, bed=args.bed, output=args.output)

# end




