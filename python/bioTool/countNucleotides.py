# coding:utf-8
def countNucleotides(seq, isDNA=True):
    if isDNA:
        print 'A:', seq.count('A')
        print 'G:', seq.count('G')
        print 'C:', seq.count('C')
        print 'T:', seq.count('T')
    else:
        print 'A:', seq.count('A')
        print 'G:', seq.count('G')
        print 'C:', seq.count('C')
        print 'U:', seq.count('U')