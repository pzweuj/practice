a = open('Broad.human.exome.b37.interval_list', 'r')
b = open('Broad.human.exome.b37.bed', 'w')

for line in a:
	if line.startswith('@'):
		continue
	else:
		lineAS = line.split('\n')[0].split('\t')
		l = [lineAS[0], lineAS[1], lineAS[2], lineAS[4], "0", lineAS[3]]
		out = '\t'.join(l)
		b.write(out + '\n')

b.close()
a.close()
