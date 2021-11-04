# coding=utf-8
# pzw
# NIPT analysis
# 20211101
# 依赖于bwa, samtools, sambamba, bedtools, R

import os
import subprocess
import argparse
import sys
import argcomplete

def main(function, option):
    now = os.path.abspath(os.path.dirname(sys.argv[0]))

    function_dict = {
        "Reference": "reference.py",
        "ZScore": "ztest.py"
    }

    print("python3 " + now + "/function/" + function_dict[function] + " " + " ".join(option))
    subprocess.call("python3 " + now + "/function/" + function_dict[function] + " " + " ".join(option), shell=True)

if __name__ == "__main__":
    now = os.path.dirname(sys.executable)
    parser = argparse.ArgumentParser(
        description="NIPT",
        prog="NIPT.py",
        usage="python3 NIPT.py [-h] <function> <function option>",
        formatter_class=argparse.RawTextHelpFormatter
    )
    parser.add_argument("-v", "--version", action="version", version="Version 0.1 20211104")
    parser.add_argument("function", choices=("Reference", "ZScore"),
        help="""
            Reference                Create NIPT Reference File
            ZScore                   Z Score
    """)
    parser.add_argument("option", nargs=argparse.REMAINDER, metavar="function option")
    argcomplete.autocomplete(parser)
    if len(sys.argv[1:]) == 0:
        parser.print_help()
        parser.exit()
    args = parser.parse_args()
    main(function=args.function, option=args.option)
