# -*- coding: utf-8 -*-
# author fengyitong 2019-01

import os
import sys
import pandas as pd
from collections import defaultdict


if len(sys.argv) != 5:
    exit('USAGE: %s infile go_level pathway out' %sys.argv[0])

info_excel = sys.argv[1]
go = sys.argv[2]
kegg = sys.argv[3]
out = sys.argv[4]

acc2go=defaultdict(set)
acc2kegg=defaultdict(set)
kegg2des = dict()

go_info = pd.read_table(go, header= 0, sep = '\t')
go_des = zip(go_info['GO'].tolist(), go_info['name'].tolist())
go2des = {x[0]:x[1] for x in go_des}
go_acc = zip(go_info['GO'].tolist(), go_info['Accs'].tolist())
for x in go_acc:
    for i in x[1].split(';'):
        acc2go[i].add(x[0])

# {acc2go[i].add(x[0]) for x in go_acc for i in x[1].split(';')}
# acc2go = {i:x[0] for x in go_acc for i in x[1].split(';')}

kegg_info = pd.read_table(kegg, header= 0, sep = '\t')
kegg_des = zip(kegg_info['pathway'].tolist(), kegg_info['pathway_name'].tolist())
kegg2des = {x[0]:x[1] for x in kegg_des}
kegg_acc = zip(kegg_info['pathway'].tolist(), kegg_info['accs_list'].tolist())
# acc2kegg = {acc2kegg[i].add(x[0]) for x in kegg_acc for i in x[1].split(';')}
for x in kegg_acc:
    for i in x[1].split(';'):
        acc2kegg[i].add(x[0])

info_pd = pd.read_excel(info_excel, header= 0, sep = '\t')
info_pd.drop(['GO','KEGG'], axis=1, inplace=True)
info_dict = info_pd.to_dict('record')
new_info_dict = list()
for dict_ in info_dict:
    acc = dict_['ID']
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

new_df = pd.DataFrame(info_dict)
new_df = new_df[info_pd.columns.tolist()+['go_info', 'kegg_info']]
new_df.to_csv(out, sep='\t',index=False,header=True)