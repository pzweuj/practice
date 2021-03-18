#!/use/bin/python3
# coding:utf-8
from selenium import webdriver
from selenium.webdriver import Edge
from bs4 import BeautifulSoup
import os
import time
import random





# HP:0000002
# HP:3000079


# hpo = "HP:0000001"
# index = 0
# maxN = 3000079

def chinahpo(hpo, dic):
    time.sleep(random.randint(5, 10))
    driver = Edge(r"C:\Program Files (x86)\Microsoft\Edge\Application\msedgedriver.exe")
    hpid = hpo.split(":")[1]
    url = "http://www.chinahpo.org/#/searchList?trigger=1&tabType=1&searchContent=HP%3A{hpid}".format(hpid=hpid)

    try:
        driver.get(url)
        strtemp = url
        print("网址：", strtemp)
    except Exception:
        print("get page error", hpo)

    time.sleep(2)
    with open("html/hp_" + hpid + ".html", "a+", encoding="utf-8") as f:
        f.write(str(driver.page_source))

    driver.close()

    dic[hpo] = {
        "url": "",
        "en_name": "",
        "cn_name": "",
        "en_def": "",
        "cn_def": ""
    }


    file = open("html/hp_" + hpid + ".html", "rb")
    html = file.read()
    soup = BeautifulSoup(html, "html.parser")
    m = soup.select("main")
    c = m[0].find_all("div", {"class": "row_list"})

    try:
        p = c[0].select("p")
        dic[hpo]["url"] = url
        dic[hpo]["en_name"] = p[0].string.split("：")[1]
        dic[hpo]["cn_name"] = p[1].string.split("：")[1]
        dic[hpo]["en_def"] = p[2].string.split("：")[1]
        dic[hpo]["cn_def"] = p[3].string.split("：")[1]
        file.close()
    except:
        print("未找到描述:", hpo)
        file.close()
        os.remove("html/hp_" + hpid + ".html")

    return dic


hpo_dict = {}
hpoFile = open("hpolist.txt", "r")
for line in hpoFile:
	hpo = line.replace("\n", "")
	hpo_dict = chinahpo(hpo, hpo_dict)
hpoFile.close()

results = open("chpo.txt", "w", encoding="utf-8")
results.write("HPID\tEn\tCn\tEn_def\tCn_def\tURL\n")

for hd in hpo_dict.keys():
    print("正在写入", hd)
    output = [hd, hpo_dict[hd]["en_name"], hpo_dict[hd]["cn_name"], hpo_dict[hd]["en_def"], hpo_dict[hd]["cn_def"], hpo_dict[hd]["url"]]
    results.write("\t".join(output) + "\n")

results.close()