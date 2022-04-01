from collections import defaultdict
import pandas as pd

org2a = dict()
with open('organisms.txt', 'r') as org_r:
    for line in org_r:
        line = line.strip().split('\t')
        org2a[line[2]] = line[0]

mi2org = dict()
with open('mirna.txt', 'r') as mir:
    for line in mir:
        line = line.strip().split('\t')
        for i in org2a:
            if i in line[4]:
                mi2org[line[2]] = i
                break
        # mi2org[line[2]] = line[4]

mi=pd.read_excel('miRNA.xls')
mi = mi.fillna('-')
mi_dict=mi.to_dict('record')
with open('miR_pre_mature.dat', 'w') as dat_w:
    dat_w.write('Species name\t3species\tpri_miR\tmature_miR_1\tmature_miR_2\n')
    for dic in mi_dict:
        org = mi2org[dic['ID']]
        orga = org2a[org]
        info = org + '\t' + orga + '\t' + dic['ID'] + '\t'
        if dic['Mature1_ID'] != '-':
            info += dic['Mature1_ID'] + '\t'
        else:
            info += '\t' + '\t'
        if dic['Mature2_ID'] != '-':
            info += dic['Mature2_ID'] + '\n'
        else:
            info += '\t' + '\n'
        dat_w.write(info)
