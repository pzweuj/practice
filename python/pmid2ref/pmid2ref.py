# pzw
# 20180212

pmidlist = open('pmid.txt', 'r')
output = open('results.txt', 'w')

def pmid2ref(pmid):
	import requests
	from bs4 import BeautifulSoup
	html = requests.get('https://www.ncbi.nlm.nih.gov/pubmed/' + str(pmid) + '/')
	soup = BeautifulSoup(html.text)
	title = soup.title.string.split('- PubMed')[0]
	info = soup.select('meta')
	for meta in info:
		if 'author' in str(meta):
			author = meta['content']
		if 'description' in str(meta):
			publish = meta['content']
			if '[' in str(publish):
				publish = publish.split('[')[0]

	results = title + ' ' + author + ' ' + publish
	return results
	
for pmid in pmidlist:
	if '\n' in pmid:
		pmid = pmid.split('\n')[0]
	# print pmid + ' done'
	results = pmid2ref(pmid).encode('utf8')
	output.write(results + '\n')

pmidlist.close()
output.close()
# print 'task done'