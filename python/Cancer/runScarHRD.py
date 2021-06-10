#!/usr/bin/python
# coding=utf-8
# pzweuj
# 20210608
# run scarHRD

import os
import sys
import argparse


def main(seg, outputDir, reference, seqz, ploidyOutputDir, chr):
    if not os.path.exists(outputDir):
        os.makedirs(outputDir)

    if chr == "TRUE" or chr == "true" or chr == "True" or chr == "T":
        chrIn = "TRUE"
    else:
        chrIn = "FALSE"

    if seqz == "TRUE" or seqz == "true" or seqz == "True" or seqz == "T":
        seqzIn = "TRUE"
    else:
        seqzIn = "FALSE"
    
    if ploidyOutputDir == "FALSE":
        cmd = """
            Rscript -e \\
                'scarHRD::scar_score("{seg}", reference="{reference}", seqz={seqzIn}, chr={chrIn})'
        """.format(seg=seg, reference=reference, seqzIn=seqzIn, chrIn=chrIn)
    else:
        cmd = """
            Rscript -e \\
                'scarHRD::scar_score("{seg}", reference="{reference}", seqz={seqzIn}, chr={chrIn}, ploidy="{ploidyOutputDir}")'
        """.format(seg=seg, reference=reference, seqzIn=seqzIn, chrIn=chrIn, ploidyOutputDir=outputDir)
    print(cmd)
    os.system(cmd)

    # results
    for i in os.listdir(os.getcwd()):
        if "_HRDresults" in i:
            resultsFile = os.getcwd() + "/" + i
            scarHRD_output = open(resultsFile, "r")
            scarHRD_output_fix = open(outputDir + "/scarHRD.results.txt", "w")
            scarHRD_output_fix.write("Sample\tHRD\tTelomeric_AI\tLST\tHRD-sum\n")
            for line in scarHRD_output:
                if not "HRD-sum" in line:
                    scarHRD_output_fix.write(line.replace('"', ''))
            scarHRD_output.close()
            scarHRD_output_fix.close()

            try:
                os.remove(resultsFile)
            except:
                pass

            try:
                os.remove(resultsFile.replace("_HRDresults.txt", "_info_seg.txt"))
            except:
                pass

if __name__ == "__main__":
    parser = argparse.ArgumentParser(
        description="scarHRD opt @pzweuj",
        prog="runScarHRD.py",
        usage="python runScarHRD.py [-h] -i <seg> -r <reference> -s <is seqz> -c <with chr> [-p <ploidy outputDir>]",
        formatter_class=argparse.RawTextHelpFormatter
    )
    parser.add_argument("-v", "--version", action="version", version="Version 0.1 20210608")
    parser.add_argument("-i", "--input", type=str, help="input file name")
    parser.add_argument("-o", "--output", type=str, help="output directory")
    parser.add_argument("-r", "--reference", type=str, help="reference genome [grch37, grch38, mouse], default: grch37", default="grch37")
    parser.add_argument("-s", "--seqz", type=str, help="if input is a seqz.gz file, default: True", default="TRUE")
    parser.add_argument("-c", "--chr", type=str, help="contains 'chr' in chromosome names, default: True", default="TRUE")
    parser.add_argument("-p", "--ploidy", type=str, help="optional, output estimated ploidy, default: False", default="FALSE")
    
    if len(sys.argv[1:]) == 0:
        parser.print_help()
        parser.exit()
    args = parser.parse_args()
    main(seg=args.input, outputDir=args.output, reference=args.reference, seqz=args.seqz, ploidyOutputDir=args.ploidy, chr=args.chr)
