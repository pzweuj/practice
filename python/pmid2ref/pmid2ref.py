# pzw
# 20180211

def pmid2ref(pmid):
	import requests
	from bs4 import BeautifulSoup
	html = requests.get('https://www.ncbi.nlm.nih.gov/pubmed/' + str(pmid) + '/')
	soup = BeautifulSoup(html.text, 'lxml')
	title = soup.title.string.split('-')[0]
	info = soup.select('meta')
	author = info[20]['content']
	publish = info[21]['content'].split('[')[0]
	results = title + ' ' + author + ' ' + publish
	return results
	
pmidlist = open('pmid.txt', 'r')
output = open('results.txt', 'w')

for pmid in pmidlist:
	if '\n' in pmid:
		pmid = pmid.split('\n')[0]
	print pmid + ' done'
	results = pmid2ref(pmid)
	output.write(results + '\n')

pmidlist.close()
output.close()
print 'task done'