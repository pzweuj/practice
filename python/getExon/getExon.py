# pzw
# 20190110

inputFile = open("ExonTest.txt", "r")
outputFile = open("exon_out.txt", "w")
top_right_num = 59595995

sequence = ""
for line in inputFile.readlines()[3:]:
	lineAfterSplit = line.split("\t")
	if len(lineAfterSplit) == 5:
		seq = "".join(lineAfterSplit[0:4])
		sequence = sequence + seq
	else:
		seq = "".join(lineAfterSplit)
		sequence = sequence + seq
sequence = sequence[::-1]

outputFile.write("seq:\n" + sequence + "\n")
outputFile.write("#Exon\tstart\tend\n")

exon = 0
count = 0
while count <= range(len(sequence) - 1):
	if sequence[count].islower():
		count += 1
	else:
		start = top_right_num + 32 - len(sequence) + count
		while count <= range(len(sequence) - 1):
			if sequence[count].isupper():
				count += 1
			else:
				exon += 1
				end = top_right_num + 31 - len(sequence) + count
				outputFile.write("exon" + str(exon) + "\t" + str(start) + "\t" + str(end) + "\n")
				break

outputFile.close()
inputFile.close()


