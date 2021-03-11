# coding=utf-8
# pzw
# 20210311

drug = open("a.txt", "r")
b = open("results.txt", "w")
drugDict = {}
for line in drug:
	if not line.startswith("Gene\t"):
		lines = line.split("\t")
		gene = lines[0]
		variant = lines[1]
		cancer = lines[2]
		significance = lines[3]
		drug = lines[4]
		level = lines[5]
		descript = lines[6].replace("\n", "")
		pmid = lines[7].split("\n")[0]

		if "+" in variant:
			variants = variant.split(" + ")
			for v in variants:
				try:
					drugDict[v].append([gene, "\t".join(v.split(" ")[1:]), drug, significance, level, cancer, descript, pmid])
				except Exception:
					drugDict[v] = []
					drugDict[v].append([gene, "\t".join(v.split(" ")[1:]), drug, significance, level, cancer, descript, pmid])
		elif "-" in variant:
			variants = variant.split(" - ")
			for v in variants:
				try:
					drugDict[v].append([gene, "fusion", drug, significance, level, cancer, descript, pmid])
				except Exception:
					drugDict[v] = []
					drugDict[v].append([gene, "fusion", drug, significance, level, cancer, descript, pmid])
		else:
			try:
				drugDict[variant].append([gene, variant.split(" ")[1], drug, significance, level, cancer, descript, pmid])
			except Exception:
				drugDict[variant] = []
				drugDict[variant].append([gene, variant.split(" ")[1], drug, significance, level, cancer, descript, pmid])

for k in drugDict.keys():
	if "mut" in k:
		continue
	if "rearrange" in k:
		continue
	if "del" in k:
		continue
	if "pos" in k:
		continue
	if "over" in k:
		continue
	if "wild" in k:
		continue
	if "amp" in k:
		continue
	if "neg" in k:
		continue
	if "loss" in k:
		continue
	if "exp" in k:
		continue
	if "hyper" in k:
		continue
	if "exon" in k:
		continue
	if "fusion" in k:
		continue

	o = []
	for i in drugDict[k]:
		o.append("[" + "|".join([i[2], i[3], i[4], i[5], i[6], i[7]]) + "]")
	output = k + "\t" + ",".join(o) + "\n"

	b.write(output)