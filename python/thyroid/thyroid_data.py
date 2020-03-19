# pzw
# 20200316

import os

def GetPanelResults(annovarFileName):
	dbCheck = open("67Panel.txt", "r")

	Panel = {}
	for line in dbCheck:
		if line.startswith("Location"):
			continue
		else:
			checkPoint = line.split("\t")
			location = checkPoint[0]
			ref = checkPoint[1]
			alt = checkPoint[2]
			ID = checkPoint[3].split("\n")[0]
			Panel[ID] = "f"

			annovarFile = open(annovarFileName, "r")
			for lineAnno in annovarFile:
				if lineAnno.startswith("Chr\tStart\tEnd"):
					continue
				else:
					c = lineAnno.split("\t")
					chrom = c[0]
					pos = c[1]
					refAnno = c[3]
					altAnno = c[4]
					DP = int(c[84].split(":")[1])
					AD = int(c[84].split(":")[2].split(",")[1])
					checkAnno = chrom + ":" + pos

					if checkAnno == location:
						if refAnno == ref and altAnno == alt:
							if AD >= 5 and DP >= 200:
								if float(AD) / float(DP) >= 0.015:
									Panel[ID] = "t"
			annovarFile.close()
	dbCheck.close()
	o = []
	for i in sorted(Panel.keys(), key=int):
		o.append(Panel[i])
	return o

def AnalysisPath(AnnovarFilePath):
	file_list = os.listdir(AnnovarFilePath)
	test = open("test_muts.arff", "w")
	test.write("@relation thyroid\n\n")
	test.write("@attribute class {Benign, Pathogenic}\n")

	dbCheck = open("67Panel.txt", "r")
	for l in dbCheck:
		if l.startswith("Loca"):
			continue
		else:
			loc = l.split("\n")[0].split("\t")[3]
			test.write("@attribute T" + loc + " { t, f}\n")
	dbCheck.close()

	test.write("\n@data\n")


	for f in file_list:
		resultsList = GetPanelResults(AnnovarFilePath + "/" + f)
		oo = "\t".join(resultsList) + "\n"
		test.write("?\t" + oo)

	test.close()

AnalysisPath("test")