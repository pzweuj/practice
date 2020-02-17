# pzw
# 20200217

import requests
from bs4 import BeautifulSoup

def ncbiNM2ensemblENST(ncbiNM):
	ncbiNM2 = ncbiNM.split(".")[0]

	html = requests.get("https://www.ncbi.nlm.nih.gov/nuccore/" + ncbiNM2)
	soup = BeautifulSoup(html.text, features="lxml")

	brief = soup.find_all(name="a", attrs={"class": "brieflinkpopperctrl"})
	for link in brief:
		if "ENST" in link.get("href"):
			ensemblENST = link.get("href").split("id/")[1]
		else:
			ensemblENST = "-"

	return ensemblENST

print ncbiNM2ensemblENST("NM_198834.1")