# coding=utf-8
# pzw
# 20200326

def Observerd_Species(inputFile, outputFile, step=1000):
	import random
	OTU = open(inputFile, "r")
	ob = open(outputFile, "w")
	ob.write("OTU_Num\tOB\tSamples\n")

	# 获得样本名
	for line in OTU:
		if line.startswith("OTU\t"):
			nameList = line.split("\n")[0].split("\t")
			nameList.pop(0)
	OTU.close()

	# 再对每个样本进行统计
	for i in range(len(nameList)):
		OTU = open(inputFile, "r")
		OTU_list = []
		for line in OTU:
			if line.startswith("OTU\t"):
				continue
			else:
				OTU_num = int(line.split("\t")[i+1])
				n = 0
				while n < OTU_num:
					OTU_list.append(line.split("\t")[0])
					n += 1
		selectNum = 0
		while selectNum < len(OTU_list):
			selectList = random.sample(OTU_list, selectNum)
			output = "\t".join([str(selectNum), str(len(set(selectList))), nameList[i]])
			ob.write(output + "\n")
			selectNum += step
		OTU.close()
	ob.close()

def shanno_index(inputFile, outputFile, step=1000):
	import random
	from collections import Counter
	import math

	OTU = open(inputFile, "r")
	sn = open(outputFile, "w")
	sn.write("OTU_Num\tshannon\tSamples\n")

	# 获得样本名
	for line in OTU:
		if line.startswith("OTU\t"):
			nameList = line.split("\n")[0].split("\t")
			nameList.pop(0)
	OTU.close()

	# 再对每个样本进行统计
	for i in range(len(nameList)):
		OTU = open(inputFile, "r")
		OTU_list = []
		for line in OTU:
			if line.startswith("OTU\t"):
				continue
			else:
				OTU_num = int(line.split("\t")[i+1])
				n = 0
				while n < OTU_num:
					OTU_list.append(line.split("\t")[0])
					n += 1
		
		## shannon index
		selectNum = 0
		while selectNum < len(OTU_list):
			selectList = random.sample(OTU_list, selectNum)
			countSelectList = Counter(selectList)
			uniqueList = list(set(selectList))
			if len(selectList) == 0:
				H = 0
			else:
				H = 0
				for u in uniqueList:
					count_u = countSelectList[u]
					pi = float(count_u) / float(len(selectList))
					# 底默认是2，一般采用2或e，如果要设置为e，则这里写math.log(pi)
					logpi = math.log(pi, 2)
					H += (logpi * pi)
				H = H * (-1)
			output = "\t".join([str(selectNum), str(H), nameList[i]])
			sn.write(output + "\n")
			selectNum += step
		OTU.close()
	sn.close()

def simpson_index(inputFile, outputFile, step=1000):
	import random
	from collections import Counter

	OTU = open(inputFile, "r")
	sp = open(outputFile, "w")
	sp.write("OTU_Num\tsimpson\tSamples\n")

	# 获得样本名
	for line in OTU:
		if line.startswith("OTU\t"):
			nameList = line.split("\n")[0].split("\t")
			nameList.pop(0)
	OTU.close()

	# 再对每个样本进行统计
	for i in range(len(nameList)):
		OTU = open(inputFile, "r")
		OTU_list = []
		for line in OTU:
			if line.startswith("OTU\t"):
				continue
			else:
				OTU_num = int(line.split("\t")[i+1])
				n = 0
				while n < OTU_num:
					OTU_list.append(line.split("\t")[0])
					n += 1
		
		## simpson index
		selectNum = 0
		while selectNum < len(OTU_list):
			selectList = random.sample(OTU_list, selectNum)
			countSelectList = Counter(selectList)
			uniqueList = list(set(selectList))
			if len(selectList) == 0:
				GS = 0
			else:
				D = 0
				for u in uniqueList:
					count_u = countSelectList[u]
					pi = float(count_u) / float(len(selectList))
					pi2 = pi ** 2
					D += pi2
				GS = 1 - D

			output = "\t".join([str(selectNum), str(GS), nameList[i]])
			sp.write(output + "\n")
			selectNum += step
		OTU.close()
	sp.close()

simpson_index("OTU.txt", "simpson_index.txt")