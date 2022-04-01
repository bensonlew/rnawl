#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""
@time    : 2019/3/28 9:09
@file    : add_info_to_proteinxlsx.py
@author  : yitong.feng
@contact: yitong.feng@majorbio.com
"""


import os
import sys
import pandas as pd
from collections import defaultdict

if not len(sys.argv) in [6,7]:
    exit('USAGE: %s infile go_level pathway diff_path description out' %sys.argv[0])

protein_excel = sys.argv[1]
go = sys.argv[2]
kegg = sys.argv[3]
diff_path = os.path.abspath(sys.argv[4])

try:
    descripton = sys.argv[5]
    out = sys.argv[6]
except:
    out = sys.argv[5]
acc2des = dict()
try:
    des_df = pd.read_csv(descripton, header=0, sep='\t')
    acc2des = {acc:des for acc,des in zip(des_df['Accession'],des_df['Description'])}
except:
    pass

out_header = list()

try:
    info_pd = pd.read_excel(protein_excel, header=0, sep='\t')
except:
    info_pd = pd.read_csv(protein_excel, header=0, sep='\t')

if acc2des:
    print('开始改变description')
    for n in range(info_pd.shape[0]):
        if info_pd.iloc[n, 3] in acc2des:
            info_pd.iloc[n, 4] = acc2des[info_pd.iloc[n, 3]]
        else:
            info_pd.iloc[n, 4] = '_'

out_header += info_pd.columns.tolist()
info_dict = info_pd.to_dict('record')

acc2diffs = defaultdict(dict)
for file in os.listdir(diff_path):
    if file.endswith('.diff.exp.xls'):
        cmp = file.split('.diff.exp.xls')[0]
        diff_file = os.path.join(diff_path, file)
        diff_info = pd.read_csv(diff_file, header=0, sep='\t')
        fc = 'FC(' + cmp.replace('_vs_', '/') + ')'
        out_header.append(fc)
        pvalue = 'Pvalue(' + cmp.replace('_vs_', '/') + ')'
        out_header.append(pvalue)
        acc2fc = zip(diff_info['Accession'].tolist(), diff_info[fc].tolist())
        for i in acc2fc:
            acc2diffs[i[0]].update({'acc':i[0],fc:i[1]})
        acc2pvalue = zip(diff_info['Accession'].tolist(), diff_info[fc].tolist())
        for i in acc2pvalue:
            acc2diffs[i[0]].update({'acc': i[0], pvalue: i[1]})

acc2go=defaultdict(set)
acc2kegg=defaultdict(set)
kegg2des = dict()

go_info = pd.read_csv(go, header= 0, sep = '\t')
go_des = zip(go_info['GO'].tolist(), go_info['name'].tolist())
go2des = {x[0]:x[1] for x in go_des}
go_acc = zip(go_info['GO'].tolist(), go_info['Accs'].tolist())
for x in go_acc:
    for i in x[1].split(';'):
        acc2go[i].add(x[0])

kegg_info = pd.read_csv(kegg, header= 0, sep = '\t')
kegg_des = zip(kegg_info['pathway'].tolist(), kegg_info['pathway_name'].tolist())
kegg2des = {x[0]:x[1] for x in kegg_des}
kegg_acc = zip(kegg_info['pathway'].tolist(), kegg_info['accs_list'].tolist())
for x in kegg_acc:
    for i in x[1].split(';'):
        acc2kegg[i].add(x[0])

diffkey = [x for x in out_header if x not in info_pd.columns.tolist()]
new_info_dict = list()
for dict_ in info_dict:
    acc = dict_['Accession']
    if acc in acc2go:
        go_info = [x + '(' + go2des[x] + ')' for x in acc2go[acc]]
        dict_.update(go_info = ';'.join(go_info))
    else:
        dict_.update(go_info='_')
    if acc in acc2kegg:
        kegg_info = [x + '(' + kegg2des[x] + ')' for x in acc2kegg[acc]]
        dict_.update(kegg_info = ';'.join(kegg_info))
    else:
        dict_.update(kegg_info='_')
    if acc in acc2diffs:
        for i in diffkey:
            if i in acc2diffs[acc]:
                dict_.update({i: str(acc2diffs[acc][i])})
            else:
                dict_.update({i: '_'})
    else:
        for i in diffkey:
            dict_.update({i: '_'})

new_df = pd.DataFrame(info_dict)
new_df = new_df[out_header + ['go_info', 'kegg_info']]
new_df.to_csv(out, sep='\t',index=False,header=True)
