a = open("a.txt", "r")
b = open("results.txt", "w")

drugDict = {}
for line in a:
	if not line.startswith("chrom"):
		lines = line.split("\t")
		chrom = lines[0]
		start = lines[1]
		end = lines[2]
		ref = lines[3]
		alt = lines[4]

		drugs = lines[33].split("&")
		drug = []
		for d in drugs:
			if "_(" in d:
				drug.append(d.split("_(")[0])
			else:
				drug.append(d)
		drug = "\\x2c".join(drug)

		signific = lines[30]
		cancer = lines[32].split("_(")[0]
		level = lines[36].split("\n")[0]
		pmid = lines[27].split("_")[0]

		checkPoint = "\t".join([chrom, start, end, ref, alt])

		try:
			drugDict[checkPoint].append([drug, signific, cancer, level, pmid])
		except Exception:
			drugDict[checkPoint] = []
			drugDict[checkPoint].append([drug, signific, cancer, level, pmid])

print(drugDict)
for i in drugDict.keys():
	o = []
	for x in drugDict[i]:
		o.append("[" + "|".join(x).replace(",", "\\x2c") + "]")
	output = i + "\t" + "\\x2c".join(o) + "\n"
	b.write(output)

