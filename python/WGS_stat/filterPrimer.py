# pzw
# 20190621

import os
import gzip
import sys

rawFastq = gzip.open(sys.argv[1], "rb")
cleandata = sys.argv[2]
cleanFastq = open(cleandata, "w")
readx = sys.argv[3]

r1adpt = "AAAGAGCATTCAAAGTGTCAAAGTAGGA"
r2adpt = "GATCGTCGGACTGTAGAACTCTGAAC"


if readx == "r1":
	adpt = r1adpt
elif readx == "r2":
	adpt = r2adpt
else:
	print "need readx"

m = 0
adaptList = []
while (m + 10) <= len(adpt):
	adaptsli = adpt[m:m+15]
	m += 1
	adaptList.append(adaptsli)


c = 0
l1 = []
l2 = []
l3 = []
l4 = []
for line in rawFastq:
	c += 1
	k = c % 4
	if k == 1:
		l1.append(line)
	elif k == 2:
		l2.append(line)
	elif k == 3:
		l3.append(line)
	elif k == 0:
		l4.append(line)
	else:
		continue

rawFastq.close()

l24 = zip(l2, l4)
l2new = []
l4new = []
for i in l24:
	cc = 0
	for ada in adaptList:
		if ada in i[0]:
			cc += 1

	if cc == 0:
		l2new.append(i[0][10:][:-10] + "\n")
		l4new.append(i[1][10:][:-10] + "\n")
	else:
		for ada in adaptList:
			if ada in i[0]:
				ii = i[0][10:].split(ada)
				cut = len(ii[0][:-10])
				l2new.append(ii[0][:-10] + ada + ii[1])
				l4new.append(i[1][10:][:cut] + i[1][10:][cut+10:])
				break

ll = zip(l1, l2new, l3, l4new)
for m in ll:
	cleanFastq.write(m[0])
	cleanFastq.write(m[1])
	cleanFastq.write(m[2])
	cleanFastq.write(m[3])





	# c += 1
	# k = c % 4
	# if k == 1:
	# 	cleanFastq.write(line)
	# elif k == 2:
	# 	line2 = line[10:]
	# 	for i in adaptList:
	# 		if i in line2:
	# 			if len(line2.split(i)[0]) >= 100:
	# 				cut = len(line2.split(i)[0][:-10])
	# 				line3 = line2.split(i)[0][:-10] + i + line2.split(i)[1]
	# 			else:
	# 				continue
	# 		else:
	# 			line3 = line2
	# 	cleanFastq.write(line3)
	# elif k == 3:
	# 	cleanFastq.write(line)
	# elif k == 0:
	# 	if cut != 0:
	# 		line4 = line[0:cut] + line[cut+10:]
	# 		cut = 0
	# 	else:
	# 		line4 = line
	# 	cleanFastq.write(line4)
	# else:
	# 	continue

cleanFastq.close()
rawFastq.close()

cmd = """
	gzip {cleandata}
""".format(cleandata=cleandata)
os.system(cmd)

