# coding:utf-8

def countPointMutation(s1, s2):
    count = 0
    for i in range(len(s1)):
        if s1[i] != s2[i]:
            count += 1
        else:
            continue
    print count