import sys
import requests
import bs4
from bs4 import BeautifulSoup
import threadpool


# 获得网页源码
def getHTMLText(url):
    try:
        r = requests.get(url, timeout=40)
        r.raise_for_status()
        r.encoding = r.apparent_encoding
        return r.text
    except:
        return ""

def print_gene_variant(links):
    url = links[-1]
    html = getHTMLText(url)
    soup = BeautifulSoup(html, "html.parser")

    des_td = []
    for td in soup.find_all("td"):
        des_td.append(td.text)

    description = "-"
    drug_res = "-"
    for des in range(len(des_td)):
        if "Gene Variant Descriptions" in des_td[des]:
            description = des_td[des + 1]
        if "Associated Drug Resistance" in des_td[des]:
            drug_res = des_td[des + 1]

    if len(soup.find_all("table", attrs={"id": "TranscriptTabTable"})) == 2:
        for tr in soup.find_all("table", attrs={"id": "TranscriptTabTable"})[1].children:
            if isinstance(tr, bs4.element.Tag):
                results = tr.text.replace(" ", "").replace("\n\n", "").strip()
                results = results.replace("\n", "\t")
                if not "GenomeBuild" in results:
                    results = results.split("\tGRCh38/hg38")[0] + "\tGRCh38/hg38"
                    results = "\t".join(links) + "\t" + description  + "\t" + drug_res + "\t" + results
                    print(results)
                    
    if len(soup.find_all("table", attrs={"id": "TranscriptTabTable"})) == 1:
        for tr in soup.find_all("table", attrs={"id": "TranscriptTabTable"})[0].children:
            if isinstance(tr, bs4.element.Tag):
                results = tr.text.replace(" ", "").replace("\n\n", "").strip()
                results = results.replace("\n", "\t")
                if not "GenomeBuild" in results:
                    results = results.split("\tGRCh38/hg38")[0] + "\tGRCh38/hg38"
                    results = "\t".join(links) + "\t" + description + "\t" + drug_res + "\t" + results
                    print(results)

    if len(soup.find_all("table", attrs={"id": "TranscriptTabTable"})) == 0:
        results = "NAN"
        results = "\t".join(links) + "\t" + description + "\t" + drug_res + "\t" + results
        print(results)

l = sys.argv[1].split("\n")[0]
links = ["test", "test", l]
print_gene_variant(links)