# pzw
# 20180212

a = open("EVTREF.fa", "r")
b = open("EVTREF_VP4.fa", "w")

seq = {}
contig = ""
c = ""
for line in a:
	if line.startswith(">"):
		c = ""
		contig = line
	else:
		c = c + line.split("\n")[0]
		seq[contig] = c

for i in seq.keys():
	seq[i] = seq[i][700: 1000]

for j in seq.keys():
	b.write(j)
	b.write(seq[j] + "\n")
