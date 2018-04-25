a = open('hg19_clinvar_20180401.vcf', 'r')
b = open('hg19_clinvar_20180401.txt', 'w')

b.write('#chrom	start	end	ref	alt	CLNSIG	CLNDN\n')

for line in a:
	if line.startswith('#'):
		continue
	else:
		lines = line.split('\t')
		chrom = lines[0]
		start = lines[1]
		ref = lines[3]
		alt = lines[4]
		end = str(int(start) + len(ref) - 1)
		clinvar = lines[7]

		clins = clinvar.split(';')

		for i in clins:
			if i.startswith('CLNDN'):
				CLNDN = i.split('=')[1]
			elif i.startswith('CLNSIG'):
				CLNSIG = i.split('=')[1]
			else:
				continue

		l = [chrom, start, end, ref, alt, CLNSIG, CLNDN]
		s = '\t'.join(l) + '\n'
		b.write(s)
b.close()
a.close()
			