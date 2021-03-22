#!/use/bin/python3
# coding:utf-8

from msedge.selenium_tools import Edge
from msedge.selenium_tools import EdgeOptions
from bs4 import BeautifulSoup
import os
import time
import random
from queue import Queue
from threading import Thread


def randomUA():
    MY_USER_AGENT = [
        "Mozilla/4.0 (compatible; MSIE 6.0; Windows NT 5.1; SV1; AcooBrowser; .NET CLR 1.1.4322; .NET CLR 2.0.50727)",
        "Mozilla/4.0 (compatible; MSIE 7.0; Windows NT 6.0; Acoo Browser; SLCC1; .NET CLR 2.0.50727; Media Center PC 5.0; .NET CLR 3.0.04506)",
        "Mozilla/4.0 (compatible; MSIE 7.0; AOL 9.5; AOLBuild 4337.35; Windows NT 5.1; .NET CLR 1.1.4322; .NET CLR 2.0.50727)",
        "Mozilla/5.0 (Windows; U; MSIE 9.0; Windows NT 9.0; en-US)",
        "Mozilla/5.0 (compatible; MSIE 9.0; Windows NT 6.1; Win64; x64; Trident/5.0; .NET CLR 3.5.30729; .NET CLR 3.0.30729; .NET CLR 2.0.50727; Media Center PC 6.0)",
        "Mozilla/5.0 (compatible; MSIE 8.0; Windows NT 6.0; Trident/4.0; WOW64; Trident/4.0; SLCC2; .NET CLR 2.0.50727; .NET CLR 3.5.30729; .NET CLR 3.0.30729; .NET CLR 1.0.3705; .NET CLR 1.1.4322)",
        "Mozilla/4.0 (compatible; MSIE 7.0b; Windows NT 5.2; .NET CLR 1.1.4322; .NET CLR 2.0.50727; InfoPath.2; .NET CLR 3.0.04506.30)",
        "Mozilla/5.0 (Windows; U; Windows NT 5.1; zh-CN) AppleWebKit/523.15 (KHTML, like Gecko, Safari/419.3) Arora/0.3 (Change: 287 c9dfb30)",
        "Mozilla/5.0 (X11; U; Linux; en-US) AppleWebKit/527+ (KHTML, like Gecko, Safari/419.3) Arora/0.6",
        "Mozilla/5.0 (Windows; U; Windows NT 5.1; en-US; rv:1.8.1.2pre) Gecko/20070215 K-Ninja/2.1.1",
        "Mozilla/5.0 (Windows; U; Windows NT 5.1; zh-CN; rv:1.9) Gecko/20080705 Firefox/3.0 Kapiko/3.0",
        "Mozilla/5.0 (X11; Linux i686; U;) Gecko/20070322 Kazehakase/0.4.5",
        "Mozilla/5.0 (X11; U; Linux i686; en-US; rv:1.9.0.8) Gecko Fedora/1.9.0.8-1.fc10 Kazehakase/0.5.6",
        "Mozilla/5.0 (Windows NT 6.1; WOW64) AppleWebKit/535.11 (KHTML, like Gecko) Chrome/17.0.963.56 Safari/535.11",
        "Mozilla/5.0 (Macintosh; Intel Mac OS X 10_7_3) AppleWebKit/535.20 (KHTML, like Gecko) Chrome/19.0.1036.7 Safari/535.20",
        "Opera/9.80 (Macintosh; Intel Mac OS X 10.6.8; U; fr) Presto/2.9.168 Version/11.52",
        "Mozilla/5.0 (Windows NT 6.1; WOW64) AppleWebKit/536.11 (KHTML, like Gecko) Chrome/20.0.1132.11 TaoBrowser/2.0 Safari/536.11",
        "Mozilla/5.0 (Windows NT 6.1; WOW64) AppleWebKit/537.1 (KHTML, like Gecko) Chrome/21.0.1180.71 Safari/537.1 LBBROWSER",
        "Mozilla/5.0 (compatible; MSIE 9.0; Windows NT 6.1; WOW64; Trident/5.0; SLCC2; .NET CLR 2.0.50727; .NET CLR 3.5.30729; .NET CLR 3.0.30729; Media Center PC 6.0; .NET4.0C; .NET4.0E; LBBROWSER)",
        "Mozilla/4.0 (compatible; MSIE 6.0; Windows NT 5.1; SV1; QQDownload 732; .NET4.0C; .NET4.0E; LBBROWSER)",
        "Mozilla/5.0 (Windows NT 6.1; WOW64) AppleWebKit/535.11 (KHTML, like Gecko) Chrome/17.0.963.84 Safari/535.11 LBBROWSER",
        "Mozilla/4.0 (compatible; MSIE 7.0; Windows NT 6.1; WOW64; Trident/5.0; SLCC2; .NET CLR 2.0.50727; .NET CLR 3.5.30729; .NET CLR 3.0.30729; Media Center PC 6.0; .NET4.0C; .NET4.0E)",
        "Mozilla/5.0 (compatible; MSIE 9.0; Windows NT 6.1; WOW64; Trident/5.0; SLCC2; .NET CLR 2.0.50727; .NET CLR 3.5.30729; .NET CLR 3.0.30729; Media Center PC 6.0; .NET4.0C; .NET4.0E; QQBrowser/7.0.3698.400)",
        "Mozilla/4.0 (compatible; MSIE 6.0; Windows NT 5.1; SV1; QQDownload 732; .NET4.0C; .NET4.0E)",
        "Mozilla/4.0 (compatible; MSIE 7.0; Windows NT 5.1; Trident/4.0; SV1; QQDownload 732; .NET4.0C; .NET4.0E; 360SE)",
        "Mozilla/4.0 (compatible; MSIE 6.0; Windows NT 5.1; SV1; QQDownload 732; .NET4.0C; .NET4.0E)",
        "Mozilla/4.0 (compatible; MSIE 7.0; Windows NT 6.1; WOW64; Trident/5.0; SLCC2; .NET CLR 2.0.50727; .NET CLR 3.5.30729; .NET CLR 3.0.30729; Media Center PC 6.0; .NET4.0C; .NET4.0E)",
        "Mozilla/5.0 (Windows NT 5.1) AppleWebKit/537.1 (KHTML, like Gecko) Chrome/21.0.1180.89 Safari/537.1",
        "Mozilla/5.0 (Windows NT 6.1; WOW64) AppleWebKit/537.1 (KHTML, like Gecko) Chrome/21.0.1180.89 Safari/537.1",
        "Mozilla/5.0 (iPad; U; CPU OS 4_2_1 like Mac OS X; zh-cn) AppleWebKit/533.17.9 (KHTML, like Gecko) Version/5.0.2 Mobile/8C148 Safari/6533.18.5",
        "Mozilla/5.0 (Windows NT 6.1; Win64; x64; rv:2.0b13pre) Gecko/20110307 Firefox/4.0b13pre",
        "Mozilla/5.0 (X11; Ubuntu; Linux x86_64; rv:16.0) Gecko/20100101 Firefox/16.0",
        "Mozilla/5.0 (Windows NT 6.1; WOW64) AppleWebKit/537.11 (KHTML, like Gecko) Chrome/23.0.1271.64 Safari/537.11",
        "Mozilla/5.0 (X11; U; Linux x86_64; zh-CN; rv:1.9.2.10) Gecko/20100922 Ubuntu/10.10 (maverick) Firefox/3.6.10",
        "Mozilla/5.0 (Windows NT 10.0; Win64; x64) AppleWebKit/537.36 (KHTML, like Gecko) Chrome/58.0.3029.110 Safari/537.36",
        "Mozilla/5.0 (Windows NT 10.0; Win64; x64) AppleWebKit/537.36 (KHTML, like Gecko) Chrome/83.0.4103.53 Safari/537.36"
    ]
    ua = random.choice(MY_USER_AGENT)
    return ua


def randomIP():
    ips = open("ip_check.txt", "r")
    ip = []
    for line in ips:
        l = line.replace("\n", "")
        ip.append(l)
    ip_random = random.choice(ip)
    return ip_random



# 爬取
def chinahpo(hpo_queue):
	
    while hpo_queue.empty() is not True:
        hpo = hpo_queue.get()

        # 如果使用IP池，则不进行随机等待
        s = random.randint(5, 10)
        print(hpo, "等待 " + str(s) + "秒")
        time.sleep(s)
        ip = randomIP()
        # ip = "socks5://127.0.0.1:1080"
        print(hpo, "使用IP " + ip)
        options = EdgeOptions()
        options.use_chromium = True
        options.add_argument("headless")
        # options.add_argument("disable-gpu")
        options.add_argument("--proxy-server={ip}".format(ip=ip))
        options.add_argument("--disable-blink-features")
        options.add_argument("--disable-blink-features=AutomationControlled")
        options.add_argument("start-maximized")
        options.add_experimental_option("excludeSwitches", ["enable-automation"])
        options.add_experimental_option("useAutomationExtension", False)
        msedge = r"C:\Program Files (x86)\Microsoft\Edge\Application\msedgedriver.exe"

        driver = Edge(options=options, executable_path=msedge)
        script = "Object.defineProperty(navigator, 'webdriver', {get: () => undefined})"
        driver.execute_script(script)
        UA = randomUA()
        # UA = "Mozilla/5.0 (Windows NT 10.0; Win64; x64) AppleWebKit/537.36 (KHTML, like Gecko) Chrome/83.0.4103.53 Safari/537.36"
        driver.execute_cdp_cmd("Network.setUserAgentOverride", {"userAgent": UA})
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
        with open("html2/hp_" + hpid + ".html", "a+", encoding="utf-8") as f:
            f.write(str(driver.page_source))

        driver.close()
        fin = open("finish.txt", "a")
        fin.write(hpo + "\n")
        fin.close()

# 解析
def analysis(hpo):
    hpid = hpo.split(":")[1]
    print(hpo)
    file = open("html/hp_" + hpid + ".html", "rb")
    html = file.read()
    soup = BeautifulSoup(html, "html.parser")
    m = soup.select("main")
    c0 = m[0].find_all("div", {"class": "el-row"})
    c = m[0].find_all("div", {"class": "row_list"})
    url = "http://www.chinahpo.org/#/searchList?trigger=1&tabType=1&searchContent=HP%3A{hpid}".format(hpid=hpid)

    output = [hpo]

    try:
        p = c[0].select("p")
        output.append(p[0].string.split("名称：")[1])
        output.append(p[1].string.split("翻译：")[1])
        output.append(p[2].string.split("定义：")[1])
        output.append(p[3].string.split("翻译：")[1])
        output.append(url)
        file.close()
    except:
        if "暂时没有该数据" in c0[1].string:
            output.append("暂时没有该数据")
        else:
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
def Grep():
    q = Queue()
    hpoFile = open("hpolist.txt", "r")
    for line in hpoFile:
        hpo = line.replace("\n", "")
        q.put(hpo)
    hpoFile.close()

    print("队列开始 %d " % q.qsize())
    for i in range(10):
        thread = Thread(target=chinahpo, args=(q, ))
        thread.start()
    q.join()




# 解析
def CreateHPO():
    hpoFileList = os.listdir("html")
    r = open("chinahpo.txt", "w", encoding="utf-8")
    r.write("HPOID\tEN\tCN\tEN_des\tCN_des\tURL\n")
    for line in hpoFileList:
        hpo = line.replace("hp_", "HP:").replace(".html", "")
        try:
            hpoString = analysis(hpo)
            r.write(hpoString + "\n")
        except:
            hpid = hpo.split(":")[1]
            filtered = open("filtered.txt", "a")
            r.write(hpo + "\t-\t-\t-\t-\thttp://www.chinahpo.org/#/searchList?trigger=1&tabType=1&searchContent=HP%3A{hpid}\n".format(hpid=hpid))
            filtered.write(hpo + "\n")
            filtered.close()
    r.close()

    filtered = open("filtered.txt", "r")
    l = []
    for line in filtered:
        if not line in l:
            l.append(line)
    filtered.close()

    filtered = open("filtered.txt", "w")
    for i in l:
        filtered.write(i)
    filtered.close()


# 重置
def remake():
    finish = open("finish.txt", "r")
    finish_list = []
    for f in finish:
        finish_list.append(f)
    finish.close()

    hpo = open("hpolist.txt", "r")
    hpolist = []
    for line in hpo:
        if not line in finish_list:
            hpolist.append(line)
    hpo.close()
    hpo = open("hpolist.txt", "w")
    for h in hpolist:
        hpo.write(h)
    hpo.close()


# 爬
Grep()

# 分析
# CreateHPO()

# 重置
# remake()