# pzw
# 20180725
# in R with dada2 use
# write.table(taxa, "taxa.txt", quote=F, sep="\t")
# write.table(seqtab.nochim, "seqtab.nochim.txt", quote=F, sep="\t")
# import pandas as pd

taxa = open('taxa.txt', 'r')
seq = open('seqtab.nochim.txt', 'r')
results = open('seq_results.txt', 'w')

dic = {}
for line in taxa:
	if line.startswith('Kingdom'):
		continue
	else:
		lineAS = line.split('\t')
		uniseq = lineAS[0]
		kingdom = lineAS[1]
		phylum = lineAS[2]
		class1 = lineAS[3]
		order = lineAS[4]
		family = lineAS[5]
		genus = lineAS[6]
		species = lineAS[7].split('\n')[0]
		l = [kingdom, phylum, class1, order, family, genus, species]
		dic[uniseq] = ';'.join(l)

newUni = []
uniChange = seq.readlines()
for i in uniChange[0].split('\n')[0].split('\t'):
	while i in dic.keys():
		i = dic[i]
		newUni.append(i)
results.write('sample\t' + '\t'.join(newUni) + '\n')

for j in uniChange[1:]:
	results.write(j)

results.close()
