# encoding:utf-8
# pzw
# 20190219
# python2.7
# v0.1
#######################
# from biomart import BiomartServer
#######################

inputdata = open("result_test_output.txt", "r")
results = open("result_anno.txt", "w")

## 查询
def searchDataset(chrom, start, end):
	from biomart import BiomartServer
	server = BiomartServer("http://asia.ensembl.org/biomart")
	server.verbose = True
	ccrigri = server.datasets["ccrigri_gene_ensembl"]
	response = ccrigri.search({
		"filters": {
			"chromosome_name": chrom,
			"start": start,
			"end": end
		},
		"attributes": [
			"ensembl_gene_id",
			"ensembl_transcript_id",
			"refseq_mrna",
			"external_gene_name",
			"start_position",
			"end_position",
			"wikigene_description"
		]
	})
	for i in response.iter_lines():
		i = i.decode("utf-8")
		return i

for line in inputdata:
	if line.startswith("hostChr"):
		continue
	lineAfterSplit = line.split("\t")
	hostChr = lineAfterSplit[0]
	hostPos1 = lineAfterSplit[1]
	hostPos2 = lineAfterSplit[2]
	searchResult = ""
	try:
		searchResult = searchDataset(hostChr, hostPos1, hostPos2)
	except Exception, ex:
		pass
	output = line.split("\n")[0] + "\t" + str(searchResult) + "\n"
	results.write(output)

results.close()
inputdata.close()