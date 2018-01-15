#!\usr\bin\python
# coding=utf-8
# pzw
# 20180115

import pandas as pd

panel = pd.read_excel('test_panel.xlsx', header=0)
amp = pd.read_excel('test_amp.xlsx', header=0)
exwriter = pd.ExcelWriter('results.xlsx')

# 生成hotspot dataframe
hotspot = pd.DataFrame(columns=['Chr', 'Start', 'End', 'UniqueID', 'Score', 'Direct', 'info', 'amp'])
hotspot['Chr'] = panel['Chr']
hotspot['Start'] = panel['Start']
hotspot['End'] = panel['End']
hotspot['UniqueID'] = panel['UniqueID']
hotspot['Score'] = '0'
hotspot['Direct'] = '+'

# 填写info列
for i in panel.index:
    if panel.loc[i, 'Alt'] == '-':
        panel.loc[i, 'Alt'] = ''
    if panel.loc[i, 'Ref'] == '-':
        panel.loc[i, 'Ref'] = ''
hotspot['info'] = 'REF=' + panel['Ref'] + ';' + 'OBS=' + panel['Alt'] + ';' + 'ANCHOR=' + panel['Anchor']

# 得到对应的引物ID
for hotspot_index in hotspot.index:
    getAmp = amp[['Amplicon_ID']][(amp['Insert_Start'] <= hotspot.at[hotspot_index, 'Start']) & (amp['Insert_Stop'] >= hotspot.at[hotspot_index, 'End']) & (amp['Chr'] == hotspot.at[hotspot_index, 'Chr'])]
    for Amp_index in getAmp.index:
        hotspot.at[hotspot_index, 'amp'] = getAmp.loc[Amp_index, 'Amplicon_ID']

# 输出hotspot到excel
hotspot.to_excel(exwriter, 'hotspot', index=False)

# 生成target的dataframe
target = pd.DataFrame(columns=['Chr', 'Start', 'End', 'amp', 'Score', 'Direct', 'info0', 'info1'])
target['Chr'] = amp['Chr']
target['Start'] = amp['Insert_Start']
target['End'] = amp['Insert_Stop']
target['amp'] = amp['Amplicon_ID']
target['Score'] = '0'
target['Direct'] = '+'
target['info0'] = '.'
target['info1'] = '.'

# 输出target到excel
target.to_excel(exwriter, 'target', index=False)
exwriter.close()

hotspot.to_csv('hotspot.bed', index=False, sep='\t', header=None)
target.to_csv('target.bed', index=False, sep='\t', header=None)