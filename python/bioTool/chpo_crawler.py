#!/use/bin/python3
# coding:utf-8

from msedge.selenium_tools import Edge
from msedge.selenium_tools import EdgeOptions
from bs4 import BeautifulSoup
import os
import time
import random

# 爬取
def chinahpo(hpo):
    time.sleep(random.randint(5, 30))
    options = EdgeOptions()
    options.use_chromium = True
    # options.add_argument("headless")
    # options.add_argument("disable-gpu")
    options.add_argument("--disable-blink-features")
    options.add_argument("--disable-blink-features=AutomationControlled")
    options.add_argument("start-maximized")
    options.add_experimental_option("excludeSwitches", ["enable-automation"])
    options.add_experimental_option("useAutomationExtension", False)
    msedge = r"C:\Program Files (x86)\Microsoft\Edge\Application\msedgedriver.exe"

    driver = Edge(options=options, executable_path=msedge)
    script = "Object.defineProperty(navigator, 'webdriver', {get: () => undefined})"
    driver.execute_script(script)
    driver.execute_cdp_cmd("Network.setUserAgentOverride", {"userAgent": "Mozilla/5.0 (Windows NT 10.0; Win64; x64) AppleWebKit/537.36 (KHTML, like Gecko) Chrome/83.0.4103.53 Safari/537.36"})
    print(driver.execute_script("return navigator.userAgent;"))

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
    fin = open("finish.txt", "a")
    fin.write(hpo + "\n")
    fin.close()

# 解析
def analysis(hpo):
    hpid = hpo.split(":")[1]
    file = open("html/hp_" + hpid + ".html", "rb")
    html = file.read()
    soup = BeautifulSoup(html, "html.parser")
    m = soup.select("main")
    c = m[0].find_all("div", {"class": "row_list"})
    url = "http://www.chinahpo.org/#/searchList?trigger=1&tabType=1&searchContent=HP%3A{hpid}".format(hpid=hpid)

    output = [hpo]

    try:
        p = c[0].select("p")
        output.append(p[0].string.split("：")[1])
        output.append(p[1].string.split("：")[1])
        output.append(p[2].string.split("：")[1])
        output.append(p[3].string.split("：")[1])
        output.append(url)
        file.close()
    except:
        output.append("未找到信息")
        output.append("-")
        output.append("-")
        output.append("-")
        output.append(url)

        # os.remove("html/hp_" + hpid + ".html")
        filtered = open("filtered.txt", "a")
        filtered.write(hpo + "\n")
        filtered.close()

    return "\t".join(output)


# 爬取
hpoFile = open("hpolist.txt", "r")
for line in hpoFile:
    hpo = line.replace("\n", "")
    chinahpo(hpo)
hpoFile.close()


# 解析
hpoFileFinish = open("finish.txt", "r")
r = open("chinahpo.txt", "w", encoding="utf-8")
r.write("HPOID\tEN\tCN\tEN_des\tCN_des\tURL\n")
for line in hpoFileFinish:
    hpo = line.replace("\n", "")
    try:
        hpoString = analysis(hpo)
        r.write(hpoString + "\n")
    except:
        pass
r.close()
hpoFileFinish.close()
