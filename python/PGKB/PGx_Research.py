# coding=utf-8
import json

inputFile = open('pagecontent-all.txt', 'r')
outputFile = open('json2.txt', 'w')

for line in inputFile:
    contents = line.split('\t')
    PGx_Research = contents[2]
    s = json.loads(PGx_Research)

    for i in range(len(s['results'])):
        l = []
        vaId = str(s['results'][i]['vaId'])
        PMID = str(s['results'][i]['lit'])
        sig = str(s['results'][i]['sig'])
        sen = '-'
        phenoCats = '-'
        genes = '-'
        races = '-'
        pValue = '-'
        cases = '-'
        chemicals = '-'
        variants = '-'
        rsid = '-'
        try:
            pValue = str(s['results'][i]['sps'][0]['pv'])
        except Exception, ex:
            pass

        try:
            genes = str(s['results'][i]['genes'][0]['geneSymbol'])
        except Exception, ex:
            pass

        try:
            sen = str(s['results'][i]['sen'])
        except Exception, ex:
            pass

        try:
            phenoCats = str(s['results'][i]['phenoCats'][0])
        except Exception, ex:
            pass

        try:
            races = str(s['results'][i]['sps'][0]['races'][0]['race'])
        except Exception, ex:
            pass

        try:
            cases = str(s['results'][i]['sps'][0]['stCase'])
        except Exception, ex:
            pass

        try:
            chemicals = ''
            for j in range(len(s['results'][i]['chemicals'])):
                chemical = s['results'][i]['chemicals'][j]['chemicalName']
                chemicals += (str(chemical) + ',')
        except Exception, ex:
            pass

        try:
            rsid = str(s['results'][i]['loc'][0]['rsid'])
        except Exception, ex:
            pass

        try:
            variants = ''
            for j in range(len(s['results'][i]['loc'])):
                variant = s['results'][i]['loc'][j]['hapName']
                variants += (str(variant) + ',')
        except Exception, ex:
            pass

        l.append(vaId)
        l.append(PMID)
        l.append(sig)
        l.append(sen)
        l.append(phenoCats)
        l.append(genes)
        l.append(races)
        l.append(pValue)
        l.append(cases)
        l.append(chemicals)
        l.append(variants)
        l.append(rsid)

        outputFile.write('\t'.join(l) + '\n')

inputFile.close()
outputFile.close()
print 'task done'