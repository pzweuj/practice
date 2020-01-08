# coding=utf-8
# pzw
# 20200108
# change ncbi origin clinvar to  annovar format

origin = open("clinvar_20200106.vcf", "r")
annovar = open("hg19_clinvar_20200106.txt", "w")

annovar.write("#Chr\tStart\tEnd\tRef\tAlt\tCLNALLELEID\tCLNDN\tCLNDISDB\tCLNREVSTAT\tCLNSIG\n")

for line in origin:
	if line.startswith("#"):
		continue
	else:
		lines = line.split("\t")
		chrom = lines[0]
		pos = lines[1]
		# ID = lines[2]
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