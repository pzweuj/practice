# coding=utf-8
# pzw
# 加密测试



import time
 
def run():
    start_time = time.time()
    n = 10000
    res = 0
    for i in range(n):
        for j in range(n):
            res += 1
    end_time = time.time()
    print("run end, use time total:", end_time - start_time)


