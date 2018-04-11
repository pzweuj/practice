# pzw
# 20180403
from bs4 import BeautifulSoup

pgkb = open('pagecontent-all.txt', 'r')
out = open('results.txt', 'w')

out.write('Vatiant	Level	Type	Gene	Phenotypes	OMB_Race	Race_Notes	Genetypes	Sen' + '\n')

for line in pgkb:
	contents = line.split('\t')
	clinical_PGx = contents[1]
	soup = BeautifulSoup(clinical_PGx, 'lxml')
	info = soup.select('.yui-gf')
	
	for i in info:
		dd = i.select('dd')
		dt = i.select('dt')
		zipped = zip(dt, dd)
		d = {
			'Level': '-',
			'Types': '-', 
			'Variant': '-', 
			'Genes': '-', 
			'Phenotypes': '-', 
			'OMB_Race': '-',
			'Race_Notes': '-', 
			'genetype': '-', 
			'sen': '-'
		}
		Level_of_evidence = zipped[0][1].string
		d['Level'] = Level_of_evidence
		
		for j in zipped:
			if j[0].string == 'Type':
				d['Types'] = j[1].string
			elif j[0].string == 'Variant':
				d['Variant'] = j[1].a.string
			elif j[0].string == 'Genes':
				d['Genes'] = j[1].a.string
			elif j[0].string == 'Phenotypes':
				d['Phenotypes'] = j[1].a.string
			elif j[0].string == 'OMB Race':
				d['OMB_Race'] = j[1].string.strip()
			elif j[0].string == 'Race Notes':
				d['Race_Notes'] = j[1].string.strip()
			else:
				continue


		span = i.select('span')
		del span[0]
		td = i.select('td')
		span_td = zip(span, td)
		for k in span_td:
			d['genetype'] = k[0].string.strip()
			d['sen'] = k[1].string

			l = [
				d['Variant'], 
				d['Level'], 
				d['Types'], 
				d['Genes'], 
				d['Phenotypes'], 
				d['OMB_Race'], 
				d['Race_Notes'], 
				d['genetype'], 
				d['sen']
			]

			output = ('\t'.join(l) + '\n').encode('utf8')
			out.write(output)


