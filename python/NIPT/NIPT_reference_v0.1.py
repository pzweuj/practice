# pzw
# 20190320

import os
import sys
import pandas
##########

path = sys.argv[1]
files = os.listdir(path)

data = {}
for file in files:
	name = file.split(".")[0]
	filex = open(path + "/" + file, "r")
	chroms = {}
	for line in filex:
		chrom = line.split("\t")[0]
		mappedR = line.split("\t")[2]
		chroms[chrom] = mappedR
	data[name] = chroms
	filex.close()

df = pandas.DataFrame(data)
df.to_csv("reference_counts.txt", sep="\t")

exit()
