a = open('TSVC_variants_IonXpress_070_071.vcf', 'r')
b = open('result.vcf', 'w')

for line in a:
	if line.startswith('#'):
		b.write(line)
	else:
		c = line.split('\t')
		chrom = c[0]
		pos = c[1]
		id = c[2]
		ref = c[3]
		alt = c[4]
		qual = c[5]
		filter = c[6]
		info = c[7]
		format = c[8]
		gt = c[9]

		id_s = id.split(';')
		alt_s = alt.split(',')

		zipped = zip(id_s, alt_s)
		for i in range(len(zipped)):
			id_new = zipped[i][0]
			alt_new = zipped[i][1]

			l = []
			l.append(chrom)
			l.append(pos)
			l.append(id_new)
			l.append(ref)
			l.append(alt_new)
			l.append(qual)
			l.append(filter)
			l.append(info)
			l.append(format)
			l.append(gt)

			d = '\t'.join(l)

			b.write(d)