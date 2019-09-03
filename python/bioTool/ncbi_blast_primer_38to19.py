# pzw
# 20180212

a = open("a.txt", "r")


seq = ""
for line in a:
	if line.startswith(">"):
		b = line.split(">")[1].split(":")[1].split("-")
		start = int(b[0])
		end = int(b[1])
	else:
		s = line.split("\n")[0]
		seq = seq + s


fp = "ACGGCTGTCCAAGGAGCTG"
rp = "GCGGATGGCGCTGAGG"

fl = len(seq.split(fp)[0])
start1 = start + fl
# print start1

rp_rev = rp[::-1].replace("A", "t").replace("T", "a").replace("C", "g").replace("G", "c").upper()
rl = len(seq.upper().split(rp_rev)[1])
end1 = end - rl

start2 = start1 + len(fp) - 1
end2 = end1 - len(rp) + 1


print "\t".join([str(start1), str(start2), str(end2), str(end1)])