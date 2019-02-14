# coding:utf-8
# pzweuj
def reverseComplement(seq, isDNA=True):
    from string import maketrans
    if isDNA:
        seq = seq.replace('U', 'T')
        transTable = maketrans('ATGCMRWSYKVHDBNatgcmrwsykvhdbn', 'TACGKYWSRMBDHVNtacgkywsrmbdhvn')
    else:
        seq = seq.replace('T', 'U')
        transTable = maketrans('AUGC', 'UACG')
    complement = seq.translate(transTable)
    reverseComp = complement[::-1]
    return  reverseComp