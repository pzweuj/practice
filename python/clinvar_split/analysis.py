# pzw
# 20180309

a = open('clinvar.txt', 'r')
b = open('results.txt', 'w')

for i in a:
	if i.startswith('#'):
		continue
	else:
		c = i.split('\n')[0]
		lineAS = c.split('\t')
		chr = lineAS[0]
		start = lineAS[1]
		end = lineAS[2]
		ref = lineAS[3]
		alt = lineAS[4]
		clinsig = lineAS[5]
		clindbn = lineAS[6]

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

		b.write('\t'.join(l) + '\n')
	
b.close()
a.close()