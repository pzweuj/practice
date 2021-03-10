from liftover import ChainFile

converter = ChainFile("/home/bioinfo/ubuntu/software/ucsc_exe/hg38ToHg19.over.chain.gz", "hg38", "hg19")

a = open("hg38.bed", "r")
b = open("hg19.bed", "w")

for line in a:
	l = line.split("\n")[0].split("\t")
	chrom = l[0]
	start = int(l[1])
	end = int(l[2])



	start_lift = converter[chrom][start][0][1]
	end_lift = converter[chrom][end][0][1]

	b.write("\t".join([chrom, str(start_lift), str(end_lift)]) + "\n")



