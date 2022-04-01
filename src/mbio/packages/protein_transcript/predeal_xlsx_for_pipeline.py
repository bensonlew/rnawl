#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""
@time    : 2019/5/21 10:23
@file    : predeal_xlsx_for_pipeline.py
@author  : yitong.feng
@contact: yitong.feng@majorbio.com
"""


import pandas as pd
import os
import argparse
import sys

reload(sys)
sys.setdefaultencoding("utf-8")

class PredealXlsx(object):
    def __init__(self, data_path, p_type, exp_flag, sample_names, db_fasta, reference, yun, out_dir, creat_des):
        self.data_path = os.path.abspath(data_path)
        if not os.path.exists(self.data_path):
            sys.exit(self.data_path+'不存在')
        self.p_type = p_type.lower()
        if self.p_type not in ['itraq', 'labelfree']:
            sys.exit('现在只能处理itraq和labelfree的原始数据')
        if not exp_flag:
            if self.p_type == 'itraq':
                exp_flag = 'Ratio:'
            else:
                exp_flag = '(Normalized):'
        self.exp_flag = exp_flag
        self.sample_names = sample_names.split(',')
        while '' in self.sample_names:
            self.sample_names.remove('')
        self.db_fasta = os.path.abspath(db_fasta)
        if not os.path.exists(self.db_fasta):
            sys.exit('搜的库不存在')
        self.reference = reference.lower()
        self.yun = yun.lower()
        if not out_dir or not os.path.exists(out_dir):
            self.out_dir = self.data_path
        else:
            self.out_dir = os.path.abspath(out_dir)
        self.creat_des = False
        if creat_des.lower() == 'yes':
            self.creat_des = True

    def get_dfs(self):
        files = os.listdir(self.data_path)
        for file in files:
            if file.endswith('.xlsx') and 'rotein' in file:
                os.rename(os.path.join(self.data_path,file),os.path.join(self.data_path,'protein.xlsx'))
            if file.endswith('.xlsx') and 'eptide' in file:
                os.rename(os.path.join(self.data_path, file), os.path.join(self.data_path, 'peptide.xlsx'))
            if file.endswith('.xlsx') and ('PSM' in file or 'psm' in file):
                os.rename(os.path.join(self.data_path, file), os.path.join(self.data_path, 'psm.xlsx'))
        protein_df = pd.read_excel(os.path.join(self.data_path,'protein.xlsx'),sep='\t')
        # 如果accession好不合规的话，需要操作一下
        for n in range(protein_df.shape[0]):
            if u'|' in protein_df.iloc[n, 3]:
                protein_df.iloc[n, 3] = protein_df.iloc[n, 3].split('|')[1]
            if u' ' in protein_df.iloc[n, 3]:
                protein_df.iloc[n, 3] = protein_df.iloc[n, 3].split(' ')[0]
        # peptide_df = pd.read_excel(os.path.join(self.data_path,'peptide.xlsx'),sep='\t')
        # psm_df = pd.read_excel(os.path.join(self.data_path,'psm.xlsx'),sep='\t')
        if self.yun == 'yes':
            peptide_df = pd.read_excel(os.path.join(self.data_path, 'peptide.xlsx'), sep='\t')
            psm_df = pd.read_excel(os.path.join(self.data_path, 'psm.xlsx'), sep='\t')
            protein_df.to_csv(os.path.join(self.out_dir,'protein.xls'), sep='\t', header=True, index=False)
            peptide_df.to_csv(os.path.join(self.out_dir,'peptide.xls'), sep='\t', header=True, index=False)
            psm_df.to_csv(os.path.join(self.out_dir,'psm.xls'), sep='\t', header=True, index=False)
        print('xlsx读取完毕')
        return protein_df

    def extract_exp(self, protein_df):
        exp_cols = list()
        # print(self.exp_flag)
        # print(protein_df.columns.tolist())

        for i in protein_df.columns:
            if self.exp_flag in i:
                exp_cols.append(i)
        exp_df = protein_df[['Accession']+exp_cols]
        exp_df = exp_df.set_index('Accession')
        # exp_df.to_csv('tmp',sep='\t',header=True,index=True)
        if self.p_type == 'labelfree':
            exp_df = exp_df.dropna(how='all', subset=exp_df.columns.tolist())
        else:
            exp_df = exp_df.dropna(how='any', subset=exp_df.columns.tolist())
            if self.reference == 'no':
                # exp_df = exp_df.insert(0,'%s_1'%self.exp_flag, pd.Series([1] * exp_df.shape[0]))
                exp_df.insert(0,'%s_1'%self.exp_flag, pd.Series([1] * len(exp_df.index.tolist()) ,index=exp_df.index))
                exp_cols.insert(0,'%s_1'%self.exp_flag)
        if len(self.sample_names) != len(exp_cols):
            print(exp_cols)
            sys.exit('样品的数量跟protein表格里面的不一样，请查看')
        col2name = {a:b for a,b in zip(exp_cols, self.sample_names)}
        exp_df = exp_df.rename(index=str, columns=col2name)
        exp_df.to_csv(os.path.join(self.out_dir,'exp.txt'), sep='\t', header=True, index=True)
        print('exp.txt生成')
        # accs = exp_df['Accession'].tolist()
        accs = exp_df.index.tolist()
        with open(os.path.join(self.out_dir,'exp.list'), 'w') as el:
            el.write('\n'.join(accs))

    def extract_fasta(self):
        cmd = 'python /mnt/ilustre/users/ting.kuang/ALL-SCRIPT/get_fasta_from_verylarge_db.py %s %s %s %s' %(os.path.join(self.out_dir,'exp.list'), self.db_fasta, 1000, os.path.join(self.out_dir,'exp.fasta'))
        os.system(cmd)
        print('exp.fasta生成')

    def extract_description(self, protein_df):
        des_df = protein_df[['Accession', 'Description']]
        des_df.to_csv(os.path.join(self.out_dir,'description.txt'), sep='\t', header=True, index=False)
        print('description.txt生成')

    def run(self):
        protein_df = self.get_dfs()
        self.extract_exp(protein_df)
        self.extract_fasta()
        if self.creat_des:
            self.extract_description(protein_df)


if __name__ == '__main__':
    # def __init__(self, data_path, p_type, exp_flag, sample_names, db_fasta, reference, yun, out_dir):
    parser = argparse.ArgumentParser(description="处理itraq和labelfree流程的原始数据")
    parser.add_argument("-data_path", type=str, required=True, help="存放的原始文件的路径")
    parser.add_argument("-p_type", type=str, default='itraq', help="项目类型，只能选itraq和labelfree两种")
    parser.add_argument("-exp_flag", type=str, default='', help="protein表格里的表达量列的共有字段，比如labelfree是(Normalized):")
    parser.add_argument("-sample_names", type=str, required=True, help=',分隔的样本名称字段，注意顺序跟protein表里的顺序一致')
    parser.add_argument("-db_fasta", type=str, required=True, help='搜库用的蛋白序列数据库')
    parser.add_argument("-reference", type=str, default='no',help='如果是itraq项目判定下有没有内参，有内参的话不会在前面多加一列1')
    parser.add_argument("-yun", type=str, default='no',help='会不会被用来跑云平台，如果是的话，会生成三个xlsx的文本文件')
    parser.add_argument("-out_dir", type=str, default='',help='输出路径，默认为空，为空的话会跟data_path一致')
    parser.add_argument("-creat_des", type=str, default='yes',help='是否生成description.txt文件')

    args = parser.parse_args()
    predeal = PredealXlsx(args.data_path, args.p_type, args.exp_flag, args.sample_names, args.db_fasta, args.reference, args.yun, args.out_dir, args.creat_des)
    predeal.run()