# coding:utf-8

def findDNAMotif(seq, key):
    con = []
    for i in range(len(seq)):
        if seq[i:i+len(key)] == key:
            con.append(str(i+1))
        else:
            continue
    if con == []:
        print 'no motif found'
    else:
        print ','.join(con)