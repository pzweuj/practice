# coding=utf-8
# pzw
# 20200108
# change ncbi origin clinvar to  annovar format

import sys
import gzip
import argparse

def main(inputs, outputs):
	origin = gzip.open(inputs, "rb")
	annovar = open(outputs, "w")

	annovar.write("#Chr\tStart\tEnd\tRef\tAlt\tCLNALLELEID\tCLNDN\tCLNDISDB\tCLNREVSTAT\tCLNSIG\n")

	for line in origin:
		if line.startswith("#"):
			continue
		else:
			lines = line.split("\t")
			chrom = lines[0]
			pos = lines[1]
			ref = lines[3]
			alt = lines[4]
			infos = lines[7].split(";")

			Chr = chrom
			Start = pos
			End = int(pos) + len(ref) - 1
			Ref = ref
			Alt = alt

			CLNALLELEID = "."
			CLNDN = "."
			CLNDISDB = "."
			CLNREVSTAT = "."
			CLNSIG = "."

			for i in infos:
				if "ALLELEID" in i:
					CLNALLELEID = i.split("=")[1]
				elif "CLNDN" in i:
					CLNDN = i.split("=")[1].replace(",", "\\x2c")
				elif "CLNDISDB" in i:
					CLNDISDB = i.split("=")[1].replace(",", "\\x2c")
				elif "CLNREVSTAT" in i:
					CLNREVSTAT = i.split("=")[1].replace(",", "\\x2c")
				elif "CLNSIG" in i:
					CLNSIG = i.split("=")[1]
				else:
					continue


			output = [Chr, Start, str(End), Ref, Alt, CLNALLELEID, CLNDN, CLNDISDB, CLNREVSTAT, CLNSIG]
			annovar.write("\t".join(output) + "\n")

	annovar.close()
	origin.close()

if __name__ == "__main__":
	parser = argparse.ArgumentParser(
		description="Clinvar TransFormar  PZW@Genephar",
		prog="clinvar_annovar_format.py",
		usage="python clinvar_annovar_format.py -i <input clinvar.gz> -o <output clinvar.txt>",
		formatter_class=argparse.RawTextHelpFormatter
	)
	parser.add_argument("-v", "--version", action="version",
		version="Version 0.2 20200109")
	parser.add_argument("-i", "--input", type=str,
		help="Input the clinvar.gz file")
	parser.add_argument("-o", "--output", type=str,
		help="output the clinvar database")

	if len(sys.argv[1:]) == 0:
		parser.print_help()
		parser.exit()
	args = parser.parse_args()
	main(inputs=args.input, outputs=args.output)