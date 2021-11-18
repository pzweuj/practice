# coding=utf-8
# pzw
# 20211107
# 解析OMIM页面

import os
from bs4 import BeautifulSoup

# 标题
def getTitle(soup, mimNum):
    title = soup.title.string.strip().split(mimNum + " - ")[1]
    return title

# 表格
def getTable(soup):
    table = soup.table
    location = phenotype = mimNumber = inheritance = mappingKey = "."
    outputString = ""
    if table:
        trs = table.find_all("tr")
        for tr in trs:
            tds = tr.find_all("td")
            tds_len = len(tds)
            if tds_len == 5:
                gene = geneNum = "."
                location = tds[0].get_text().strip() if tds[0].get_text().strip() != "" else "."
                phenotype = tds[1].get_text().strip() if tds[1].get_text().strip() != "" else "."
                mimNumber = tds[2].get_text().strip() if tds[2].get_text().strip() != "" else "."
                inheritance = tds[3].get_text().strip() if tds[3].get_text().strip() != "" else "."
                mappingKey = tds[4].get_text().strip() if tds[4].get_text().strip() != "" else "."
            elif tds_len == 4:
                location = gene = geneNum = "."
                phenotype = tds[0].get_text().strip() if tds[0].get_text().strip() != "" else "."
                mimNumber = tds[1].get_text().strip() if tds[1].get_text().strip() != "" else "."
                inheritance = tds[2].get_text().strip() if tds[2].get_text().strip() != "" else "."
                mappingKey = tds[3].get_text().strip() if tds[3].get_text().strip() != "" else "."
            elif tds_len == 7:
                location = tds[0].get_text().strip() if tds[0].get_text().strip() != "" else "."
                phenotype = tds[1].get_text().strip() if tds[1].get_text().strip() != "" else "."
                mimNumber = tds[2].get_text().strip() if tds[2].get_text().strip() != "" else "."
                inheritance = tds[3].get_text().strip() if tds[3].get_text().strip() != "" else "."
                mappingKey = tds[4].get_text().strip() if tds[4].get_text().strip() != "" else "."
                gene = tds[5].get_text().strip() if tds[5].get_text().strip() != "" else "."
                geneNum = tds[6].get_text().strip() if tds[6].get_text().strip() != "" else "."
            else:
                continue
            outputString = outputString + "[" + "|".join([location, phenotype, mimNumber, inheritance, mappingKey, gene, geneNum]) +"]##"
    return outputString.rstrip("##") if outputString != "" else "."

# clinical synopsis
def getClinicalFold(soup, filter=False):
    clinicalSynopsisFoldList = soup.select("#clinicalSynopsisFold")
    childStringDiv = "."
    if len(clinicalSynopsisFoldList) != 0:
        clinicalDivList = clinicalSynopsisFoldList[0].select(".small")[0]
        childDiv = clinicalDivList.find_all("div", recursive=False)
        childStringDiv = ""
        for i in range(len(childDiv) - 1):
            childDivList = childDiv[i].get_text().replace("\n", "").strip().split(" - ")
            headerDiv = childDivList[0].strip()
            if filter:
                if headerDiv == "INHERITANCE":
                    continue
                if headerDiv == "MOLECULAR BASIS":
                    continue
                if headerDiv == "Inheritance":
                    continue
            stringDiv = ""
            for i in range(len(childDivList)):
                if i != 0:
                    if filter:
                        stringDiv = stringDiv + childDivList[i].strip().split(" [")[0] + ", "
                    else:
                        stringDiv = stringDiv + childDivList[i].strip() + ", "
            if filter:
                stringDiv = stringDiv.rstrip(", ")
            else:
                stringDiv = headerDiv + ": " + stringDiv.rstrip(", ")
            childStringDiv = childStringDiv + stringDiv + "; "
        childStringDiv = childStringDiv.rstrip("; ")
    return childStringDiv

# description
def getDescription(soup):
    descriptionFoldList = soup.select("#descriptionFold")
    descriptionFold = "." if len(descriptionFoldList) == 0 else descriptionFoldList[0].get_text().strip().replace("\n", " \\ ")
    return descriptionFold

# combine
def main(mimGene, inputDir, outputFile):
    mimGenes = open(mimGene, "r", encoding="utf-8")
    mimGenesDict = {}
    for line in mimGenes:
        if line.startswith("#"):
            continue
        elif line.startswith("ID"):
            continue
        else:
            lines = line.split("\t")
            id = lines[0]
            gene = lines[3]
            mimGenesDict[id] = gene
    mimGenes.close()

    outputs = open(outputFile, "w", encoding="utf-8")
    outputs.write("#MimNum\tGene\tTitle\tTable\tDescription\tClinical\tFilter\tURL\n")
    for i in os.listdir(inputDir):
        if ".html" in i:
            print("正在解析：", i)
            mimNum = i.split(".html")[0]
            gene = mimGenesDict[mimNum]
            url = "https://omim.org/entry/{}".format(mimNum)
            html = open(inputDir + "/{}.html".format(mimNum), "r", encoding="utf-8")
            soup = BeautifulSoup(html, "lxml")
            html.close()
            title = getTitle(soup, mimNum)
            table = getTable(soup)
            if gene == "":
                if table == ".":
                    gene = "."
                else:
                    genes = table.split("##")
                    gene_list = []
                    for g in genes:
                        gg = g.split("|")[5]
                        if not gg == ".":
                            gene_list.append(gg)
                    if len(gene_list) == 0:
                        gene = "."
                    else:
                        gene = ",".join(gene_list)

            clinical = getClinicalFold(soup)
            clinical_filt = getClinicalFold(soup, True)
            des = getDescription(soup)
            output = "\t".join([mimNum, gene, title, table, des, clinical, clinical_filt, url])
            outputs.write(output + "\n")
    outputs.close()

main("mim2gene.txt", "entry", "test.txt")
