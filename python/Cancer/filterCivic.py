# coding=utf-8
# pzw
# 20210310


civic = open("nightly-civic_accepted_and_submitted.vcf", "r")
civic_r = open("results.txt", "w")

for line in civic:
	if not line.startswith("#"):
		lineAfterSplit = line.split("\t")
		chrom = "chr" + lineAfterSplit[0]
		start = lineAfterSplit[1]
		ref = lineAfterSplit[3]
		alt = lineAfterSplit[4]

		if len(ref) > len(alt):
			t = "del"
			end = str(int(start) + len(ref) - 1)
		elif len(ref) < len(alt):
			t = "ins"
			end = start
		else:
			t = "normal"
			end = start

		Info = lineAfterSplit[7]
		Infos = Info.split(";")

		Infos_dict = {}
		for i in Infos:
			k = i.split("=")[0]
			v = i.split("=")[1]
			Infos_dict[k] = v


		gene = Infos_dict["GN"]
		variant = Infos_dict["VT"]
		CSQ = Infos_dict["CSQ"].split("\n")[0]
		print(CSQ)
		if "," in CSQ:
			c = CSQ.split(",")
			for i in c:
				ii = i.replace("|", "\t")
				o = "\t".join([chrom, start, end, ref, alt, gene, variant]) + "\t" + ii + "\n"
				civic_r.write(o)
		else:
			CSQ = CSQ.replace("|", "\t")
			o = "\t".join([chrom, start, end, ref, alt, gene, variant]) + "\t" + CSQ + "\n"
			civic_r.write(o)

