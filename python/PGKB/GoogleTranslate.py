# coding=utf-8
import urllib
import urllib2

text = 'Hello'
values={'hl':'zh-CN','ie':'UTF-8','text':text,'langpair':"en|zh-CN"}
url='http://translate.google.cn/translate_t'
data = urllib.urlencode(values)
req = urllib2.Request(url, data)
req.add_header('User-Agent', "Mozilla/4.0 (compatible; MSIE 6.0; Windows NT 5.1; SV1; .NET CLR 2.0.50727)")
response = urllib2.urlopen(req)

print response.read()