# pzw
# 20200319

import os

file_list = os.listdir("test")
outputFile = open("test.matrix.txt", "w")
changeList = []

for files in file_list:
	f = open("test/" + files, "r")
	for line in f:
		if line.startswith("Chr\tStart"):
			continue
		else:
			l = line.split("\t")
			chrom = l[0]
			start = l[1]
			end = l[2]
			exonicOrNot = l[5]
			gene = l[6]
			nonsynon = l[8]
			changeInfo = l[9]

			if gene == "TERT":
				print changeInfo
			else:
				if exonicOrNot != "exonic":
					continue
				elif "nonsynonymous" not in nonsynon:
					continue
				else:
					cc = changeInfo.split(":")
					transcript = cc[1]
					exon = cc[2]
					baseChange = cc[3]
					aachange = cc[4]

					if "," in aachange:
						aachange = aachange.split(",")[0]

					output = "\t".join([gene, chrom + ":" + start + "-" + end, transcript, exon, baseChange, aachange])
					changeList.append(output)
	f.close()

changeList = sorted(list(set(changeList)))

line0 = "changeID"
line1 = "Gene"
line2 = "Location"
line3 = "transcript"
line4 = "exon"
line5 = "baseChange"
line6 = "ID/aaChange"

n = 0
for i in changeList:
	n += 1
	line0 = line0 + "\t" + str(n)
	ii = i.split("\t")
	line1 = line1 + "\t" + ii[0]
	line2 = line2 + "\t" + ii[1]
	line3 = line3 + "\t" + ii[2]
	line4 = line4 + "\t" + ii[3]
	line5 = line5 + "\t" + ii[4]
	line6 = line6 + "\t" + ii[5]


outputFile.write(line0 + "\n")
outputFile.write(line1 + "\n")
outputFile.write(line2 + "\n")
outputFile.write(line3 + "\n")
outputFile.write(line4 + "\n")
outputFile.write(line5 + "\n")
outputFile.write(line6 + "\n")

for files in file_list:
	f = open("test/" + files, "r")

	Results_f = []
	for insertInfo in changeList:
		insertX = ""
		fx = open("test/" + files, "r")
		for line in fx:
			if line.startswith("Chr\tStart"):
				continue
			else:
				l = line.split("\t")
				chrom = l[0]
				start = l[1]
				end = l[2]
				nonsynon = l[8]
				changeInfo = l[9]
				otherInfo = l[84]
				DP = otherInfo.split(":")[1]
				AF = otherInfo.split(":")[2].split(",")[1]

				checkPoint = chrom + ":" + start + "-" + end

				if checkPoint == insertInfo.split("\t")[1]:
					if int(AF) > 10 and int(DP) > 100:
						if float(AF) / float(DP) > 0.015:
							MAF =  "%.2f" % (float(AF) / float(DP))
							insertX = MAF + "(" + AF + "/" + DP + ")"
		fx.close()
		Results_f.append(insertX)
	outputFile.write(files.split(".")[0] + "\t" + "\t".join(Results_f) + "\n")

	f.close()
outputFile.close()