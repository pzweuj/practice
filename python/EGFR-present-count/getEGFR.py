# pzw
# 20180109

filename = open('filename.txt', 'r')
EGFRconclusion = open('existEGFR.txt', 'w')


def existEGFR(inputfile, Hotspot=True):
    import pandas as pd
    df = pd.read_csv(inputfile, header=0, sep='\t')
    if Hotspot == True:
    	df2 = df.loc[(df['Allele Source'] == 'Hotspot') & (df['Allele Call'].str.contains('H'))]
    else:
    	df2 = df.loc[df['Allele Call'].str.contains('H')]
    l = []
    for i in df2['Gene ID']:
        l.append(i)
    if 'EGFR' in l:
        conclue = inputfile + '-True'
    else:
        conclue = inputfile + '-False'
    del df
    del df2
    return conclue

for name in filename:
    openFilePath = './inputdata/' + name.split('\n')[0]
    temp = existEGFR(openFilePath)
    EGFRconclusion.write(temp + '\n')

EGFRconclusion.close()
filename.close()