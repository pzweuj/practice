# pzw
# 20190927

# ref: https://www.ncbi.nlm.nih.gov/pmc/articles/PMC4951315/bin/oncotarget-07-18631-s001.pdf


def OrBaseRiskCalculation(resultsDict, database):
	db = open(database, "r")
	dbDict = {}
	for line in db:
		if line.startswith("#"):
			continue
		else:
			lines = line.split("\n")[0].split("\t")
			chrom = lines[2]
			pos = lines[3]
			alt = lines[5]
			OR = float(lines[10])
			frequence = float(lines[9])
			poskey = chrom + ":" + pos
			w = (frequence ** 2) * (OR ** 2) + 2 * frequence * (1 - frequence) * OR + ((1 - frequence) ** 2)
			dbDict[poskey] = [alt, OR, w]
	db.close()
	score_total = 1
	score_normal = 1
	for i in resultsDict.keys():
		n = resultsDict[i].count(dbDict[i][0])
		score = (float(dbDict[i][1]) ** n) / dbDict[i][2]
		score_normal = score_normal * (1 / dbDict[i][2])
		score_total = score_total * score

	risk_power = score_total / score_normal

	return score_total #, score_normal, risk_power

res = {
	"chr6:31540784": "CC",
	"chr5:52347369": "CT",
	"chr9:22125503": "GC",
	"chr17:61565891": "--",
}


print OrBaseRiskCalculation(res, "CVDRiskDB.txt")
