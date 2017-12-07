# -*- coding:utf-8 -*-
from bs4 import BeautifulSoup

inputFile = open('pagecontent-all.txt', 'r')
outputFile = open('html.txt', 'w')

for line in inputFile:
    contents = line.split('\t')
    clinical_PGx = contents[1]
    soup = BeautifulSoup(clinical_PGx, 'lxml')

    for i in range(len(soup.select('.yui-gf'))):
        temp = soup.select('.yui-gf')[i]
        dd = temp.select('dd')
        dt = temp.select('dt')
        d = {'Level': '-', 'Types': '-', 'Variant': '-', 'Genes': '-', 'Phenotypes': '-', 'OMB_Race': '-','Race_Notes': '-'}
        zipped = zip(dt, dd)
        Level_of_Evidence = zipped[0][1].string
        d['Level'] = Level_of_Evidence
        for j in range(len(zipped)):
            if j == 0:
                continue
            stan = zipped[j][0].em.string

            if stan == 'Type':
                d['Types'] = zipped[j][1].string
            elif stan == 'Variant':
                variants = zipped[j][1].select('a')
                m = []
                for n in range(len(variants)):
                    vars = variants[n].string
                    m.append(vars)
                d['Variant'] = ','.join(m)
            elif stan == 'Genes':
                d['Genes'] = zipped[j][1].a.string
            elif stan == 'OMB Race':
                d['OMB_Race'] = zipped[j][1].string.strip()
            elif stan == 'Race Notes':
                d['Race_Notes'] = zipped[j][1].string.strip()
            elif stan == 'Phenotypes':
                d['Phenotypes'] = zipped[j][1].a.string
            else:
                continue

            for k in range(len(temp.select('span'))):
                l = []
                if k == 0:
                    continue
                else:
                    genetype = temp.select('span')[k].string.strip()
                    sen = temp.select('td')[k - 1].string

                l.append(d['Variant'])
                l.append(d['Level'])
                l.append(d['Types'])
                l.append(d['Genes'])
                l.append(d['Phenotypes'])
                l.append(d['OMB_Race'])
                l.append(d['Race_Notes'])
                l.append(genetype)
                l.append(sen)

                outputFile.write(('\t'.join(l) + '\n').encode('utf-8'))

inputFile.close()
outputFile.close()
print 'task done'