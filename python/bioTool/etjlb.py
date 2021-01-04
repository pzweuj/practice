# 20210104
"""
  儿童俱乐部
X        儿
------------
部部部部部部

求每个字代表的数字
"""

# for etjlb in range(10000, 100000):
#     e = etjlb // 10000
#     t = etjlb // 1000 - (e * 10)
#     j = etjlb // 100 - (e * 100) - (t * 10)
#     l = etjlb // 10 - (e * 1000) - (t * 100) - (j * 10)
#     b = etjlb - (e * 10000) - (t * 1000) - (j * 100) - (l * 10)
    
#     bbbbbb = b + b * 10 + b * 100 + b * 1000 + b * 10000 + b * 100000
#     if (etjlb * e) == bbbbbb:
#         print(etjlb)

for etjlb in range(10000, 100000):
    etjlb_list = []
    for i in str(etjlb):
        etjlb_list.append(i)
    e = int(etjlb_list[0])
    b = etjlb_list[4]

    bbbbbb = int(b + b + b + b + b + b)
    b = int(b)
    if (etjlb * e) == bbbbbb:
        print(etjlb)