# pzw
# 20190806

from bs4 import BeautifulSoup
import json

soup = BeautifulSoup(open("show.php"), "html.parser")

for i in soup.find_all("script"):
	ageList = []
	carrierList = []
	nonCarrierList = []
	if "dash" in str(i.string):
		a = i.string.split("FusionCharts(")[1].split(");")[0]

		j = json.loads(a)

		sig = j["dataSource"]["chart"]["caption"]
		print sig

		type_age = j["dataSource"]["categories"][0]["category"]
		carrier_risk = j["dataSource"]["dataset"][0]["data"]
		nonCarrier_risk = j["dataSource"]["dataset"][1]["data"]

		for j in type_age:
			ageList.append(j["label"])

		for k in carrier_risk:
			carrierList.append(k["value"])

		for m in nonCarrier_risk:
			nonCarrierList.append(m["value"])

		Zipped = zip(ageList, carrierList, nonCarrierList)

		for n in Zipped:
			print n[0], n[1], n[2]