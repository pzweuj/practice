# coding=utf-8
# 20211107

import requests
import json
import random
import time

# ["mimNumber", "entryType", "geneSymbol"]
def readMim2Gene(filename):
    fileList = []
    with open(filename, "r") as f:
        for line in f.readlines():
            if not line.startswith("#"):
                strs = line.split("\t")
                mimNumber = strs[0]
                mimEntryType = strs[1]
                geneSymbol = "."
                if strs[3] != "":
                    geneSymbol = strs[3]
                if not "moved" in mimEntryType:
                    tempList = [mimNumber, mimEntryType, geneSymbol]
                    fileList.append(tempList)
    return fileList

# IP pool
# https://raw.fastgit.org/fate0/proxylist/master/proxy.list
def ip_pool(listFile):
    p = open(listFile, "r")
    proxy_list = []
    for line in p:
        j = json.loads(line)
        host = j["host"]
        port = j["port"]
        anonymity = j["anonymity"]

        if anonymity == "high_anonymous":
            output = host + ":" + str(port)
            proxy_list.append(output)
    p.close()
    return proxy_list

# spider
def spider(mimNumber, ip_use):
    UA = [
        "Mozilla/5.0 (compatible; bingbot/2.0; +http://www.bing.com/bingbot.htm)",
        "Mozilla/5.0 AppleWebKit/537.36 (KHTML, like Gecko; compatible; bingbot/2.0; +http://www.bing.com/bingbot.htm) Chrome/96.0.4664.33 Safari/537.36 Edg/95.0.1020.40"
    ]
    url = "https://omim.org/entry/{}".format(mimNumber)
    proxy = {"https": "https://" + ip_use, "http": "http://" + ip_use}
    ua_use = random.choice(UA)
    header = {"User-Agent": ua_use}
    s = requests.session()
    if ip_use != "":
        s.proxies = proxy
    s.headers = header
    s.keep_alive = False
    
    try:
        html = s.get(url)    
        with open("entry/" + mimNumber + ".html", "w", encoding="utf-8") as f:
            f.write(html.text)
        ip_use_log = open("ip_use.txt", "a")
        ip_use_log.write(mimNumber + "\t" + ip_use + "\n")
        ip_use_log.close()
    except:
        print(mimNumber, "无法连接，无法获取")

# main
def running(mim2gene, ipFile):
    if ipFile != "":
        ip = ip_pool(ipFile)
        ip_use = random.choice(ip)
    else:
        ip_use = ""
    i = open("ip_use.txt", "r")
    idList = []
    for line in i:
        ids = line.replace("\n", "").split("\t")[0]
        idList.append(ids)
    mimList = readMim2Gene(mim2gene)
    i.close()
    for m in mimList:
        mimNumber = m[0]
        if not mimNumber in idList:
            print("正在爬取", mimNumber, ip_use)
            spider(mimNumber, ip_use)
            sleepTime = random.randint(4, 10)
            print("等待", sleepTime, "秒")
            time.sleep(sleepTime)

running("mim2gene.txt", "proxy.list")
