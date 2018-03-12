# pzw
# 20180309

clinvar = open('clinvar.txt', 'r')
clinvarSplit = open('clinvarSplit.txt', 'w')

for i in clinvar:
	if i.startswith('#'):
		continue
	else:
		iSplitN = i.split('\n')[0]
		line = iSplitN.split('\t')
		chr = line[0]
		start = line[1]
		end = line[2]
		ref = line[3]
		alt = line[4]
		clinsig = line[5]
		clindbn = line[6]

	clinsigList = clinsig.split('|')
	clindbnList = clindbn.split('|')

	zipped = zip(clinsigList, clindbnList)
	for ii in zipped:
		clsig = ii[0]
		cldbn = ii[1]

		# print clsig

		l = []
		l.append(chr)
		l.append(start)
		l.append(end)
		l.append(ref)
		l.append(alt)
		l.append(clsig)
		l.append(cldbn)

		clinvarSplit.write('\t'.join(l) + '\n')
	
clinvarSplit.close()
clinvar.close()