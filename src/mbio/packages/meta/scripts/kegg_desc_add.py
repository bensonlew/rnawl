# -*- coding:utf-8 -*-

import pandas as pd
import sys
import json
import copy
import argparse
import re
import os
a = argparse.ArgumentParser()
a.add_argument('-k',help='xls of KEGG K. KO ID in head')
a.add_argument('-L2', help='xls of kegg level2. CategoryL2 in head')
a.add_argument('-L3', help='xls of kegg level3. CategoryL3 in head')
a.add_argument('-k_database' , help='database contain: K id,Level3, Level2, Level1,PATHWAY')
a.add_argument('-k_json', help='json file,contain K information')
args = a.parse_args()

k_file = args.k   #'predictions_ko.xls'
k2_file = args.L2  #'predictions_ko.L2.xls'
k3_file = args.L3  #'predictions_ko.L3.xls'
database_file = args.k_database  #'/mnt/ilustre/users/sanger-dev/app/bioinfo/meta/16sFundb/database/ko.category.txt'  #sys.argv[5]
k_json_file = args.k_json   #'/mnt/ilustre/users/sanger-dev/app/bioinfo/meta/16sFundb/database/kegg_ko.detail.json'

with open(k_json_file) as k_json_handle:
    aa = k_json_handle.readline()
    # k_json = json.load(k_json_handle)
    k_json = eval(aa)

#add K desc
k_data_ori = pd.read_table(k_file,sep='\t', header=0)

k_data = copy.deepcopy(k_data_ori)

def _get_k_des(x):
    if x in k_json:
        return  k_json[str(x)]['ko_desc']
    else:
        return '--'

new_columns = list(k_data.columns)
k_data['K_desc'] = k_data['KO ID'].apply(lambda x: _get_k_des(x))
os.rename(k_file, k_file+'_old')
new_columns.insert(1, 'K_desc')
k_data.to_csv(k_file,sep='\t',index=False, columns=new_columns)

## produce module

module_desc = {}
replace_pat = re.compile('\[PATH:[map\d\s]*\]')
def _get_module(x):
    if x not in k_json:
        return ''
    if k_json[x]['module_id']:
        modules = ';'.join(k_json[x]['module_id'])
        for id, m in enumerate(k_json[x]['module_id']):
            if m not in module_desc:
                print(x)
                print(modules)
                print(k_json[x]['module_category'])
                tmp_desc =  k_json[x]['module_category'][id][-1]
                tmp_desc = re.sub(replace_pat,'',tmp_desc)
                module_desc[m] = tmp_desc
    else:
        modules = ''
    return modules

k_data_ori['Module'] = k_data_ori['KO ID'].apply(lambda x: _get_module(str(x)))
k_data_ori = k_data_ori.drop('Module',axis=1).join(k_data_ori['Module'].str.split(';',expand=True).stack().reset_index(level=1,drop=True).rename('Module'))
k_data_ori=k_data_ori.reset_index()
k_data_ori.drop('KO ID',axis=1)

module_stat = k_data_ori.groupby('Module').sum()
module_stat = module_stat.reset_index()
module_stat = module_stat.drop(0,axis=0)  #去除总计

def _get_module_desc(x):
    if x in module_desc:
        return module_desc[x]
    else:
        return '--'



new_columns = list(module_stat.columns)
module_stat['description'] = module_stat['Module'].apply(lambda x: _get_module_desc(x))
module_stat = module_stat[module_stat['Module'] != "-"]

new_columns.insert(1,'description')
module_file = 'predictions_module.xls'
module_stat.to_csv(module_file, sep='\t', index=False, columns=new_columns)


ref_data = pd.read_table(database_file, sep='\t', header=0)
#level2 add level1
k2_data = pd.read_table(k2_file, sep='\t', header=0)

ex_data_2 = ref_data[['Level2', 'Level1']]
ex_data_2.drop_duplicates(inplace=True)
k2_new_data = pd.merge(k2_data,ex_data_2,how='left',left_on='CategoryL2', right_on='Level2')

k2_new_data.fillna("-", inplace=True)
new_columns = list(k2_data.columns)
new_columns.insert(1,'Level1')
os.rename(k2_file, k2_file+"_old")
k2_new_data.to_csv(k2_file,sep='\t',index=False,columns=new_columns)


#level3 and level2 , level1
k3_data = pd.read_table(k3_file, sep='\t', header=0)

ex_data_3 = ref_data[['Level3', 'Level2', 'Level1','PATHWAY']]
ex_data_3['PATHWAY'] = ex_data_3['PATHWAY'].apply(lambda x: str(x).replace('PATH:',''))
ex_data_3.drop_duplicates(inplace=True)
k3_new_data = pd.merge(k3_data,ex_data_3,how='left',left_on='CategoryL3', right_on='Level3')
k3_new_data['PATHWAY'] = k3_new_data['PATHWAY'].apply(lambda x: str(x).replace('nan','-'))
k3_new_data.sort_values(by="PATHWAY", ascending=False, inplace=True)
k3_new_data.drop_duplicates(subset=['CategoryL3'], keep='first', inplace=True)
k3_new_data.fillna("-", inplace=True)

new_columns = list(k3_data.columns)
new_columns.insert(1,'PATHWAY')
new_columns.insert(2,'Level2')
new_columns.insert(3,'Level1')

os.rename(k3_file, k3_file+"_old")
k3_new_data.to_csv(k3_file,sep='\t',index=False,columns=new_columns)

