import re

gDNA = open('loc_results_hg38.txt', 'r')
results = open('results.txt', 'w')

for line in gDNA:
	chr = line.split(':')[0]
	loc0 = line.split(':')[1]
	loc1 = loc0.split('.')[1]
	
	if loc1.__contains__('ins'):
		start = loc1.split('_')[0]
		end_t = loc1.split('_')[1]
		end = end_t.split('ins')[0]
		ref = '-'
		alt = end_t.split('ins')[1]

		# print chr, start, end, ref, alt
	elif loc1.__contains__('del'):
		start = loc1.split('del')[0]
		alt = '-\n'
		ref = loc1.split('del')[1].split('\n')[0]
		end = str(int(start) + len(ref))

	else:
		loc2 = loc1.split('>')
		alt = loc2[1]
		l = re.findall(r'\d+|\D', loc2[0])
		start = l[0]
		end = start
		ref = l[1]

	result = []
	result.append(chr)
	result.append(start)
	result.append(end)
	result.append(ref)
	result.append(alt)

	results.write('\t'.join(result))

gDNA.close()
results.close()
