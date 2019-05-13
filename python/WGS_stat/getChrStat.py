# coding=utf-8
# pzw
# 20190513

sample = 
outputfile = open("results.txt", "w")

chromList = ["chrX", "chrY", "chrM"]
i = 1
while i <= 22:
	chrom = "chr" + str(i)
	i += 1
	chromList.append(chrom)

def getChromCoverage(inputfile, chrom):
	inputFile = open(inputfile, "r")
	length = 0
	covLength = 0
	counts = 0
	for line in inputFile:
		lineAS = line.split("\t")
		chromosome = lineAS[0]
		start = lineAS[1]
		end = lineAS[2]
		count = lineAS[3].split("\n")[0]
		if chromosome == chrom:
			length_tmp = int(end) - int(start)
			length += length_tmp
			counts += (int(count) * length_tmp)
			if count != "0":
				covLength_tmp = length_tmp
				covLength += covLength_tmp

	coverage = "%.2f" % ((float(covLength) / float(length)) * 100) + "%"
	depth_cov = "%.2f" % (float(counts) / float(covLength))
	depth_all = "%.2f" % (float(counts) / float(length))
	inputFile.close()

	return [coverage, depth_cov, depth_all]

# print getChromCoverage(sample, "chrM")

coverage_list = ["chromosome"]
depth_cov_list = ["覆盖区域的平均深度"]
depth_all_list = ["WGS平均深度"]

for m in chromList:
	chromi = m.split("\n")[0]
	print chromi, getChromCoverage(sample, chromi)
	coverage_list.append(getChromCoverage(sample, chromi)[0])
	depth_cov_list.append(getChromCoverage(sample, chromi)[1])
	depth_all_list.append(getChromCoverage(sample, chromi)[2])

outputfile.write("\t".join(chromList) + "\n")
outputfile.write("\t".join(coverage_list) + "\n")
outputfile.write("\t".join(depth_cov_list) + "\n")
outputfile.write("\t".join(depth_all_list) + "\n")

outputfile.close()



