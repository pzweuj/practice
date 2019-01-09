# pzw
# 20190109
# cruzdb/sqlalchemy/MySQL-python


snplist = open("snplist.txt", "r")
outputFile = open("snpinfo_results.txt", "w")

def getSnpInfo(x):
	from cruzdb import Genome
	hg19 = Genome(db="hg19")
	snp151 = hg19.snp151
	info = snp151.filter_by(name=x).first()
	return info

outputFile.write("#chr	start	end	id\n")

for line in snplist:
	rsid = line.split("\n")[0]
	result = getSnpInfo(rsid)
	outputFile.write(str(result) + "\n")
	print rsid

print "task done"

outputFile.close()
snplist.close()
