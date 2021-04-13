# coding=utf-8
# pzw
# 20210413

from bs4 import BeautifulSoup
import os

results = open("DECIPHER.txt", "w", encoding="utf-8")

n = 1
while n <= 55:
    htmlFile = "html/decipher_page" + str(n) + ".html"
    soup = BeautifulSoup(open(htmlFile), "lxml")

    table = soup.find_all("table")[0]
    tbody = table.select("tbody")[0]
    tr = tbody.find_all("tr")

    for t in tr:
        td = t.find_all("td")
        
        name_col = td[0].find_all("div")
        gene = name_col[0].string
        description = name_col[1].string
        
        loca_col = td[1].div.div.find_all("span")
        chrom = "chr" + loca_col[0].string
        start = loca_col[2].string
        end = loca_col[4].string

        pLI = td[2].span.string
        LOEUF = td[3].span.string
        sHet = td[4].span.string
        HI_precent = td[5].span.string
        
        OMIM_col = td[6]
        if len(OMIM_col.find_all("div")) == 0:
            OMIM = "-"
        else:
            if not "OMIM" in OMIM_col.div.a.string:
                OMIM = "-"
            else:
                try:
                    OMIM = OMIM_col.div.a["href"]
                except:
                    OMIM = "-"

        DDG2P_col = td[7]
        if not DDG2P_col.string == "-":
            DDG2P = DDG2P_col.a.div.string.strip()
        else:
            DDG2P = DDG2P_col.string

        ClinGen_col = td[8]
        if not ClinGen_col.string == "-":
            ClinGen_col = ClinGen_col.a.div.find_all("div")
            ClinGen_list = []
            for c in ClinGen_col:
                ClinGen_list.append(c.text.strip().replace("\n", "").replace(" ", "").replace("Tr", ";Tr"))
            ClinGen = ";".join(ClinGen_list)
        else:
            ClinGen = ClinGen_col.string

        Open = td[9].string
        
        Link_col = td[10].find_all("li")
        link_list = []
        for l in Link_col:
            link_list.append(l.a["href"])
        links = ",".join(link_list)
        
        results.write("\t".join([gene, description, chrom, start, end, pLI, LOEUF, sHet, HI_precent, OMIM, DDG2P, ClinGen, Open, links]) + "\n")
    n += 1

results.close()