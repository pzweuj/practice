import sys
HAinput = open(sys.argv[1], "r")
NAinput = open(sys.argv[2], "r")

outputdata = open(sys.argv[3], "w")

HA = []
dictHA = {}
for line2 in HAinput:
	if line2.startswith(">"):
		HA.append(line2.split("\n")[0].split("> ")[1])

for key2 in HA:
	dictHA[key2] = dictHA.get(key2, 0) + 1

print dictHA


NA = []
dictNA = {}
for line1 in NAinput:
	if line1.startswith(">"):
		NA.append(line1.split("\n")[0].split("> ")[1])

for key1 in NA:
	dictNA[key1] = dictNA.get(key1, 0) + 1

print dictNA

sumHA = 0
for i in dictHA:
	sumHA += dictHA[i]

sumNA = 0
for j in dictNA:
	sumNA += dictNA[j]

HAa = max(dictHA, key=dictHA.get)
HAb = "%.2f" % (float(dictHA[HAa]) / float(sumHA) * 100.00)
HAout = "HA subtype: " + HAa + ", percentage: " + HAb + "%"
print HAout

NAa = max(dictNA, key=dictNA.get)
NAb = "%.2f" % (float(dictNA[NAa]) / float(sumNA) * 100.00)
NAout = "NA subtype: " + NAa + ", percentage: " + NAb + "%"
print NAout

outputdata.write(str(dictHA) + "\n")
outputdata.write(str(dictNA) + "\n")
outputdata.write(HAout + "\n")
outputdata.write(NAout + "\n")

outputdata.close()
HAinput.close()
NAinput.close()


exit()