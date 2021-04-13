# coding=utf-8
# pzw
# 20210413
# https://www.deciphergenomics.org/genes


from msedge.selenium_tools import Edge
from msedge.selenium_tools import EdgeOptions
import time
import random

options = EdgeOptions()
options.use_chromium = True
options.add_argument("headless")
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

url = "https://www.deciphergenomics.org/genes"

driver.get(url)
print("网址:", url)
# 等待加载
time.sleep(40)

# 定位下拉选择框并选择100
driver.find_element_by_xpath('//*[@id="content"]/div/div/div[2]/div/div/div[2]/div/div[1]/div/label/select/option[@value="100"]').click()
time.sleep(10)

# 保存第一页
n = 1
f = open("html/decipher_page" + str(n) + ".html", "wb")
f.write(driver.page_source.encode("gbk", "ignore"))
print("已完成第" + str(n) + "页")
f.close()

# 总共有55页
n += 1
while n <= 55:
    t = random.randint(5, 10)
    print("等待", t, "秒")
    time.sleep(t)
    driver.find_element_by_xpath("//a[text()='{page}']".format(page=str(n))).click()
    time.sleep(5)
    f = open("html/decipher_page" + str(n) + ".html", "wb")
    f.write(driver.page_source.encode("gbk", "ignore"))
    print("已完成第" + str(n) + "页")
    f.close()
    n += 1


