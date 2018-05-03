#!/usr/bin/env python
# -*- coding:utf-8 -*-

import os
import sys
import gzip
import pandas as pd

__doc__ = '''
USAGE:
    python script <in:clinvar.vcf.gz> <out:clinvar.tab.xls> <out:annovar.annot.db>

clinvar.vcf.gz:
    archive_2.0
    ftp://ftp.ncbi.nlm.nih.gov/pub/clinvar/vcf_GRCh37/clinvar.vcf.gz

clinvar.tab.xls:
    full vcf infor

annovar.annot.db:
    #Chr	Start	End	Ref	Alt	CLINSIG	CLNDBN	CLNACC	CLNDSDB	CLNDSDBID
    1	949523	949523	C	T	Pathogenic	Immunodeficiency_38_with_basal_ganglia_calcification	RCV000162196.3	MedGen:OMIM	CN221808:616126
    1	949608	949608	G	A	Benign	not_specified	RCV000455759.1	MedGen	CN169374
    1	949696	949696	-	G	Pathogenic	Immunodeficiency_38_with_basal_ganglia_calcification	RCV000148989.5	MedGen:OMIM	CN221808:616126
'''

try:
    _, clinvar, fulltab, annovar = sys.argv
except:
    print(__doc__)
    sys.exit(1)

def parse_variant(line):
    CHROM, POS, ID, REF, ALT, QUAL, FILTER, INFO = line.split('\t')
    POS = int(POS)
    dat = {}
    dat['CHROM'] = CHROM
    dat['ID'] = ID
    if len(REF)==len(ALT)==1 and REF!=ALT:    # SNV
        dat['START'], dat['END'] = POS, POS
        dat['REF'], dat['ALT'] = REF, ALT
    elif len(REF)==1 and len(ALT)!=1 and ALT.startswith(REF):    # Insert
        dat['START'], dat['END'] = POS, POS
        dat['REF'], dat['ALT'] = '-', ALT[1:]
    elif len(REF)!=1 and len(ALT)==1 and REF.startswith(ALT):    # Delete
        dat['START'], dat['END'] = POS+1, POS+len(REF)-1
        dat['REF'], dat['ALT'] = REF[1:], '-'
    else:
        #print(line)
        dat['START'], dat['END'] = POS, POS+len(REF)-1
        dat['REF'], dat['ALT'] = REF, ALT
    INFO = [i.split('=') for i in INFO.split(';')]
    INFO = {i[0]:i[1] for i in INFO}
    dat.update(INFO)
    dat['CLNSIG'] = dat.get('CLNSIG', 'not provided').replace('/', '|')
    dat['CLNACC'] = ''
    try:
        # example  CLNDISDB=MedGen:C3808739,OMIM:615120|MedGen:CN169374
        tmp = [i.split(',') for i in dat['CLNDISDB'].split('|')]
        db, ids = [], []
        for i in tmp:
            tmp2 = []
            for j in i:
                tmp2.append(['.', '.'] if j=='.' else j.split(':'))
            db.append(','.join([j[0] for j in tmp2]))
            ids.append(','.join([j[1] for j in tmp2]))
        dat['CLNDSDB'] = '|'.join(db)
        dat['CLNDSDBID'] = '|'.join(ids)
    except KeyError:
        pass
    return pd.Series(dat)


fopen = gzip.open if clinvar.endswith('gz') else open
with fopen(clinvar) as f:
    vcf = []
    for i in f:
        if i.startswith('#') or not i.strip():
            continue
        vcf.append(parse_variant(i.strip()))

dat = pd.concat(vcf, axis=1, join='outer')
dat = dat.T
#dat = dat.fillna('--')
title = ['CHROM', 'START', 'END', 'REF', 'ALT', 'ID']
header = sorted(list(set(dat.columns)-set(title)))
dat = dat[title+header]
dat.to_csv(fulltab, sep='\t', index=False)

# header = ['CHROM', 'START', 'END', 'REF', 'ALT', 'CLNSIG', 'CLNDN', 'CLNACC', 'CLNDSDB', 'CLNDSDBID']
header = ['CHROM', 'START', 'END', 'REF', 'ALT', 'CLNSIG', 'CLNDN', 'CLNDISDB']
dat = dat[header]
header = ['#Chr', 'Start', 'End', 'Ref', 'Alt', 'CLINSIG', 'CLNDBN', 'CLNDISDB']
dat.columns = header
dat.to_csv(annovar, sep='\t', index=False)
