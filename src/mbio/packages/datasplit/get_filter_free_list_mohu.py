#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""
@time    : 2019/11/26 12:43
@file    : get_filter_free_list.py
@author  : yitong.feng
@contact: yitong.feng@majorbio.com
"""

import pandas as pd
import numpy as np
import os
from collections import defaultdict
import argparse




def get_unique_index(df, inter, use='Description'):
    unique_dict = dict()
    perfer_dict = dict()
    des_ = use
    for ind in df.index:
        des = str(df.loc[ind, des_])
        com_id = str(df.loc[ind, 'Compound ID'])
        # 还得忽略一下大小写
        des = des.lower()
        if inter:
            if len(des) - len(inter) > 3:
                continue
        score = df.loc[ind, 'Fragmentation Score']
        try:
            if com_id.startswith('HMDB') or score < 30:
                pass
        except Exception:
            print(com_id, df, ind)
        if com_id.startswith('HMDB') or score < 30:
            perfer = 0
        else:
            perfer = 1
        if not des in perfer_dict:
            perfer_dict[des] = perfer
        ind_ = unique_dict.get(des)
        if not ind_:
            unique_dict[des] = ind
        else:
            perfer_ = perfer_dict.get(des)
            if perfer_ and not perfer:
                # 如果本次是HMDB的就不用比了
                continue
            if perfer and not perfer_:
                # 如果上次是HMDB而上一次不是则直接改掉
                unique_dict[des] = ind
                perfer_dict[des] = perfer
                continue
            score_ = df.loc[ind_, 'Fragmentation Score']
            if score_ < score:
                unique_dict[des] = ind
    unique_dict = {v: k for k, v in unique_dict.items()}
    unique_index = list(unique_dict.keys())
    unique_index.sort()
    if inter == 'chlorogenic acid':
        print(inter, unique_index)
    return unique_index



def get_filter_free_list(pos_i, neg_i, inter_f, out_f, out_f_n):
    # raw_dir = os.getcwd()
    # os.chdir(os.path.dirname(pos_i))
    with open(pos_i, 'rb') as pir:
        pos_i_df = pd.read_csv(pir, sep=',', quoting=0, quotechar='"', encoding='utf-8', usecols=['Compound', 'Compound ID', 'Description', 'Fragmentation Score']).fillna(0)
    with open(neg_i, 'rb') as nir:
        neg_i_df = pd.read_csv(nir, sep=',', quoting=0, quotechar='"', encoding='utf-8', usecols=['Compound', 'Compound ID', 'Description', 'Fragmentation Score']).fillna(0)


    dfs = list()
    with open(inter_f) as inter_i, open(out_f, 'wb') as ow, open(out_f+'_except.list', 'w') as ew:
        for inter in inter_i:
            inter = inter.strip().split('\t')[0]
            if not inter:
                continue
            # # pos_inter = pos_i_df[pos_i_df['Description'].replace(np.nan,'0').astype(str).str.contains(inter[1:])]
            # pos_inter = pos_i_df[pos_i_df['Description'].replace(np.nan,'0').astype(str).str.contains(inter,regex=False)]
            # # neg_inter = neg_i_df[neg_i_df['Description'].replace(np.nan,'0').astype(str).str.contains(inter[1:])]
            # neg_inter = neg_i_df[neg_i_df['Description'].replace(np.nan,'0').astype(str).str.contains(inter,regex=False)]
            try:
                pos_inter = pos_i_df[pos_i_df['Description'].replace(np.nan,'0').astype(str).str.contains(inter[1:], case=False)]
            except Exception:
                pos_inter = pos_i_df[pos_i_df['Description'].replace(np.nan,'0').astype(str).str.contains(inter,regex=False)]
            try:
                neg_inter = neg_i_df[neg_i_df['Description'].replace(np.nan,'0').astype(str).str.contains(inter[1:], case=False)]
            except Exception:
                neg_inter = neg_i_df[neg_i_df['Description'].replace(np.nan,'0').astype(str).str.contains(inter,regex=False)]
            #根据封一统20211026更新的脚本更改
            if not pos_inter.shape[0] and not neg_inter.shape[0]:
                print(inter+'在本项目中没有鉴定到')
                ew.write(inter + '\n')
                continue
            inter_dict = defaultdict(set)
            if pos_inter.shape[0]:
                pos_inter.index = list(range(0, pos_inter.shape[0]))
                # print(pos_inter.index)
                pos_inter = pos_inter.loc[get_unique_index(pos_inter, inter), :]
                for n, row in pos_inter.iterrows():
                    inter_dict[row['Description']].add(row['Compound'])
                dfs.append(pos_inter)
            if neg_inter.shape[0]:
                neg_inter.index = list(range(0, neg_inter.shape[0]))
                # print(neg_inter.index)
                neg_inter = neg_inter.loc[get_unique_index(neg_inter, inter), :]
                for n, row in neg_inter.iterrows():
                    inter_dict[row['Description']].add(row['Compound'])
                dfs.append(neg_inter)
            for des, com in inter_dict.items():
                o = des + '\n' + '\n'.join(list(com)) + '\n'
                ow.write(o.encode('utf-8'))

    df = pd.concat(dfs)
    df.index = list(range(0, df.shape[0]))
    df.to_excel(out_f+'.xlsx', encoding='utf-8', index=False)
    df = df.loc[get_unique_index(df, None, 'Compound'), :]
    inter_dict = defaultdict(set)
    for n, row in df.iterrows():
        inter_dict[row['Description']].add(row['Compound'])
    # with open('filter_free_no_dup.list', 'wb') as fw:
    with open(out_f_n, 'wb') as fw:
        for des, com in inter_dict.items():
            o = des + '\n' + '\n'.join(list(com)) + '\n'
            fw.write(o.encode('utf-8'))
    # os.chdir(raw_dir)


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description="")
    parser.add_argument("-pos", "--pos_i", dest="pos_i", required=True, type=str, help="")
    parser.add_argument("-neg", "--neg_i",dest="neg_i", required=True, type=str, help="")
    parser.add_argument("-inter", "--inter_f",dest="inter_f", required=True, type=str, help="")
    parser.add_argument("-out", "--out_f",dest="out_f", required=True, type=str, help="")
    parser.add_argument("-outn", "--out_f_n",dest="out_f_n", required=True, type=str, help="")
    args = parser.parse_args()
    get_filter_free_list(args.pos_i, args.neg_i, args.inter_f, args.out_f, args.out_f_n)
