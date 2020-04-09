# coding=utf-8
# pzw
# 20200409

import os
import sys
import argparse

def getAbsPath():
	now = os.path.abspath(os.path.dirname(sys.argv[0]))
	return now

def getLoca(loca, base):
	loc = loca.split(":")[1]
	location_start = str(int(loc) - base)
	location_end = str(int(loc) + base)
	chrom = loca.split(":")[0].split("chr")[1]
	location = chrom + ":" + location_start + "-" + location_end
	return loc, location

def main(bamfile, output, loca, row, base):
	locs = getLoca(loca, base)
	loc = locs[0]
	location = locs[1]
	bam2raster = getAbsPath() + "/function/bam2raster.jar"
	reference = "/home/genephar/databases/hg19/hg19.fa"
	cmd = """
		java -jar {bam2raster} -o {output} \\
			-r {location} -R {reference} \\
			{bamfile} --limit {row} --highlight {loc}
	""".format(bam2raster=bam2raster, output=output, location=location, reference=reference, bamfile=bamfile, loc=loc, row=row)
	os.system(cmd)
	print("Task done!")

# main("H04938Y4.final.bam", "test.png", "chr1:11174400")

if __name__ == "__main__":
	parser = argparse.ArgumentParser(
		description="Bamfile Snapshot @PZW",
		prog="snapshot.py",
		usage="python3 snapshot.py -b <bamfile> -o <png image> -l <location>",
		formatter_class=argparse.RawTextHelpFormatter
	)
	parser.add_argument("-v", "--version", action="version",
		version="Version 0.1 20200409")
	parser.add_argument("-b", "--bam", type=str,
		help="input the hg19 aligned bam file")
	parser.add_argument("-o", "--output", type=str,
		help="output png format image")
	parser.add_argument("-l", "--location", type=str,
		help="exp: '-l chr1:11174400'")
	parser.add_argument("-row", "--row", type=str,
		help="显示行数， 默认50", default="50")
	parser.add_argument("-base", "--base", type=int,
		help="前后碱基数，默认50", default=50)

	if len(sys.argv[1:]) == 0:
		parser.print_help()
		parser.exit()
	args = parser.parse_args()
	main(bamfile=args.bam, output=args.output, loca=args.location, row=args.row, base=args.base)