# pzw
# 20171205
# coding:utf-8
from bs4 import BeautifulSoup

inputFile = open('AD.html', 'r')
outputFile = open('result.txt', 'w')

soup = BeautifulSoup(inputFile, 'lxml')
gene = soup.select('p')

genes = []
for i in range(len(gene)):
    if str(gene[i]).__contains__('span class'):
        genes.append(gene[i])

genes.pop()
geneall = []
for m in range(len(genes)):
    l = str(genes[m]).split('\n')
    line = ''.join(l)
    lines = ''.join(line.split())
    genelist = lines.split('<b><i>')
    del genelist[0]
    genelist[-1] = genelist[-1].replace('</p>', '')
    for n in range(len(genelist)):
        geneall.append(genelist[n])

for x in geneall:
    genename = x.split('</i></b>:')[0]
    aachange = x.split('</i></b>:')[1]
    aachanges = aachange.split('span')
    del aachanges[0]
    del aachanges[-1]
    aachangess = ''.join(aachanges).split('>)(<')
    for item in aachangess:
        genedes = genename + '\t' + item.split('">')[0].split('"')[1] + '\t' + item.split('">')[1].split('<')[0]
        outputFile.write(genedes + '\n')

outputFile.close()
inputFile.close()
print 'task done'