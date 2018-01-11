# pzw
# 20180111

inputfile = open('../Thyroid_mutation.json', 'r')
outputfile = open('results.txt', 'w')

for line in inputfile:
	if line.startswith('  "genomic_dna_change":'):
		a = line.split('  "genomic_dna_change": "')[1]
		b = a.split('",')[0]
		outputfile.write(b + '\n')

inputfile.close()
outputfile.close()