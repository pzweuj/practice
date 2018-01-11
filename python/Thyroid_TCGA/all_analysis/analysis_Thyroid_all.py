# pzw
# 20180111
import json

inputfile = open('../Thyroid_mutation.json', 'r')
outputfile = open('results.txt', 'w')

analysisJson = json.load(inputfile)

for i in analysisJson:
    consequence = i['consequence']
    for j in consequence:
        outList = []
        transcript = j['transcript']
        outDict = {'genomic_dna_change': '-', 'mutation_subtype': '-', 'gene': '-', 'aa_change': '-', 'impact': '-'}
        outDict['genomic_dna_change'] = i['genomic_dna_change']
        outDict['mutation_subtype'] = i['mutation_subtype']
        outDict['gene'] = transcript['gene']['symbol']

        if transcript['aa_change'] != None:
            outDict['aa_change'] = transcript['aa_change']

        if transcript.has_key('annotation'):
            outDict['impact'] = transcript['annotation']['vep_impact']

        outList.append(outDict['genomic_dna_change'])
        outList.append(outDict['mutation_subtype'])
        outList.append(outDict['gene'])
        outList.append(outDict['aa_change'])
        outList.append(outDict['impact'])

        outStr = '\t'.join(outList) + '\n'

        outputfile.write(outStr)

outputfile.close()
inputfile.close()
print('task done')