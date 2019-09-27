# pzw
# 20190927

results = {
	"chr1:11854476": "AA",
	"chr1:11856378": "CC",
	"chr1:169510348": "GG",
	"chr1:169519049": "GA"
}

def calculateRisk(inputData, database):
	db = open(database, "r")
	i = 0
	snpRisk = 1
	for line in db:
		i += 1
		if i == 1:
			print "kit= ", line.split("# ")[1]
		elif i == 2:
			r = float(line.split("# ")[1].split("\n")[0])
		elif i == 3:
			continue
		elif i >= 4:
			info = line.split("\n")[0].split("\t")
			chrom = info[1]
			pos = info[2]
			alt = info[4]
			OR = info[5]
			p = info[6]
			matchKey = chrom + ":" + pos
			n = inputData[matchKey].count(alt)
			singleOR = (float(OR) ** n) / ((1 + (float(OR) - 1) * float(p)) ** 2)
			snpRisk = snpRisk * singleOR

	db.close()
	risk = r * snpRisk
	return risk, snpRisk

a = calculateRisk(results, "CVDRiskDB.txt")
print a