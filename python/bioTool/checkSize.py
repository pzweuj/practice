# coding=utf-8
# pzw
# 20210323

import os
import shutil

def getDocSize(path):
	try:
		size = os.path.getsize(path)
		return size
	except Exception as err:
		print(err)


path = "html2"
for i in os.listdir(path):
	f = path + "/" + i
	s = getDocSize(f)
	
	if s >= 9000:
		if s<= 15000:
			shutil.copy(f, "html3/" + i)