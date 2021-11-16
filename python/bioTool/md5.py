# coding=utf-8

import os
import hashlib
import sys

def get_file_md5(fname):
    m = hashlib.md5()
    with open(fname, "rb") as fobj:
        while True:
            data = fobj.read(4096)
            if not data:
                break
            m.update(data)
    return m.hexdigest()

def printMD5(path):
    for i in os.listdir(path):
        if os.path.isfile(path + "/" + i):
            f = path + "/" + i
            m = get_file_md5(f)
            print(m + "\t" + i)

p = sys.argv[1]
printMD5(p)
