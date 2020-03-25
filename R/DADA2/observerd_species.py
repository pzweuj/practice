# pzw

import random
OTU = open("count.txt", "r")

OTU_list = []
for line in OTU:
	if line.startswith("Representative_Sequence"):
		continue
	else:
		ls = line.split("\t")
		OTU_id = ls[0]
		K1 = ls[7]
		K1_num = int(K1)

		n = 0
		while n < K1_num:
			OTU_list.append(OTU_id)
			n += 1


selectNum = 0
while selectNum < len(OTU_list):
	selectList = random.sample(OTU_list, selectNum)
	selectNum += 1000
	print(len(set(selectList)))
