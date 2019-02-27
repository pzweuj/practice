# pzw
# 20190215
import sys

#############################
plasmidOrigin = open(sys.argv[1], "r")
plasmidFasta = open(sys.argv[2], "w")
plasmidName = sys.argv[3]
#############################

def getSeqLine(x, y=10):
	seq_tmp = x[y:].split(" ")
	seq = "".join(seq_tmp).upper()
	return seq

plasmidFasta.write(">" + plasmidName + "\n")
for line in plasmidOrigin:
	if line == "":
		continue
	seqs = getSeqLine(line)
	plasmidFasta.write(seqs)

plasmidFasta.close()
plasmidOrigin.close()