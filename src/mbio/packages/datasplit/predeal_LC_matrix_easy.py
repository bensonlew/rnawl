#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""
@time    : 2019/8/28 14:20
@file    : predeal_LC_matrix.py
@author  : yitong.feng
@contact: yitong.feng@majorbio.com
"""


import argparse
import os
import sys
import pickle
import pandas as pd
import copy
import chardet
from io import StringIO
from adjust_multi_CID import DealMultiCid # 20210527加上显红写的处理多个cid的功能


class PredealLcMatrix(object):
    def __init__(self, pos_i, pos_q, neg_i, neg_q, min_row, max_row, min_score, max_score, sure_score,\
        max_len, group_table, chname, cutoff_percent, filter_free, inter_compound, database_dir, output,keep_hmdb='yes',kegg_amino_acid='no'):
        self.dict_al2al = {
            u'α':u'alpha',
            u'β':u'beta',
            u'γ':u'gamma',
            u'δ':u'delta',
            u'ε':u'epsilon',
            u'ζ':u'zeta',
            u'η':u'eta',
            u'θ':u'theta',
            u'ι':u'iota',
            u'κ':u'kappa',
            u'λ':u'lambda',
            u'μ':u'mu',
            u'ν':u'nu',
            u'ξ':u'xi',
            u'ο':u'omicron',
            u'π':u'pi',
            u'ρ':u'rho',
            u'σ':u'sigma',
            u'τ':u'tau',
            u'υ':u'upsilon',
            u'φ':u'phi',
            u'χ':u'chi',
            u'ψ':u'psi',
            u'ω':u'omega',
            u'&plusmn;':u'±',
            u'&':u'',
            u';':u'',
            u'′':u"'",
            u'<WBR>':u'',
            # u'+/-': u'±',
            u'±': u'+/-',
            u'Δ': u'delta',
            u'- ': u'-',
        }

        # 用于将氨基酸右旋改为左旋
        self.amino_acid = {
            'C00037': 'Glycine',
            'C00041': 'Alanine',
            'C00183': 'Valine',
            'C00123': 'Leucine',
            'C00407': 'Isoleucine',
            'C00049': 'Aspartic acid',
            'C00152': 'Asparagine',
            'C00025': 'Glutamic acid',
            'C00064': 'Glutamine',
            'C00065': 'Serine',
            'C00188': 'Threonine',
            'C00073': 'Methionine',
            'C00097': 'Cysteine',
            'C00047': 'Lysine',
            'C00062': 'Arginine',
            'C00135': 'Histidine',
            'C00148': 'Proline',
            'C00079': 'Phenylalanine',
            'C00082': 'Tyrosine',
            'C00078': 'Tryptophan',
        }
        self.amino_acid = {v: k for k, v in self.amino_acid.items()}
        # 传入四个表格和最大行数与最小得分
        self.pos_i = os.path.abspath(pos_i)
        self.pos_q = os.path.abspath(pos_q)
        self.neg_i = os.path.abspath(neg_i)
        self.neg_q = os.path.abspath(neg_q)
        self.min_row = min_row
        self.max_row = max_row
        self.min_score = min_score
        self.max_score = max_score
        self.sure_score = sure_score
        self.max_len = max_len
        self.keep_hmdb=keep_hmdb
        self.kegg_amino_acid = kegg_amino_acid
        # 参考库
        meta_qc_pk = os.path.join(database_dir, "meta_c.pk")
        meta_qc_pk_ = os.path.join(database_dir, "CID2NAME.20210409.v98.0.txt")
        meta_hmdb_pk = os.path.join(database_dir, "meta_hmdb.pk")
        meta_hmdb_pk_ = os.path.join(database_dir, "hmdb_metabolites.detail.xls")
        extron_meta = os.path.join(database_dir, "extron_meta.txt")
        # anno_info_pk = os.path.join(database_dir, "anno_info.pk")# 更改为 -new
        anno_info_pk = os.path.join(database_dir, "anno_info_new.pk") #20220124更新 --薛钦文，封一统
        # basic_info_v10 = os.path.join(database_dir, "MJDB_selfbuild.basic_info_v10.xls")
        basic_info_v10 = os.path.join(database_dir, "MJDB.ID_Name.v12.20211130.compatible.v10.xls")
        smpdb_file = os.path.join(database_dir, "smpdb_metabolites.20180914.xls")#20220124更新 --薛钦文，封一统
        self.cid2name = os.path.join(database_dir, "CID2NAME.20210409.v98.0.txt")
        self.cid2map = os.path.join(database_dir, "Compound_Pathway.list")
        self.cid2detail = os.path.join(database_dir, "KEGG_compound.20210526_addHMDB.xls")
        self.kegg = os.path.join(database_dir, "CID2NAME.20210409.v98.0.txt")
        self.hmdb = os.path.join(database_dir, "hmdb_metabolites.detail.xls")
        self.lipids = os.path.join(database_dir, "LIPIDMAPS.detail.xls")
        self.output = output
        # self.epa = os.path.join(database_dir, "EPA.detail.xls")

        # 所有氨基酸的3字简写对1字简写的字典
        self.peps = {'Val': 'V', 'Xaa': 'X', 'Cys': 'C', 'Asp': 'D', 'Asx': 'B', 'Phe': 'F', 'Met': 'M', 'Pyl': 'O', 'Leu': 'L', 'Asn': 'N', 'Tyr': 'Y', 'Ile': 'I', 'Gln': 'Q', 'Xle': 'J', 'Thr': 'T', 'Glx': 'Z', 'Gly': 'G', 'His': 'H', 'Trp': 'W', 'Glu': 'E', 'Ser': 'S', 'Lys': 'K', 'Pro': 'P', 'Ala': 'A', 'Sec': 'U', 'Arg': 'R'}
        # 读入有kegg C号的代谢物
        # meta_c = 'Y:\\蛋白代谢事业部\\蛋白与代谢生信分析部\\06实验室研发项目\\sx_do_in_win\\代谢数据库\\meta_c.pk'
        meta_c = meta_qc_pk
        if not os.path.exists(meta_c):
            self.meta_c_dict = dict()
            # with open('Y:\\蛋白代谢事业部\\蛋白与代谢生信分析部\\06实验室研发项目\\sx_do_in_win\\代谢数据库\\CID2NAME.20210409.v98.0.txt')\
            with open(meta_qc_pk_)\
                as knr:
                for line in knr:
                    line = line.strip().split('\t')
                    if not line or len(line) < 2:
                        continue
                    for c in line[1].split(';'):
                        self.meta_c_dict[self.get_stardard(c.strip())] = 1
                        self.meta_c_dict[c.strip()] = 1
            with open(meta_c, 'wb') as mw:
                pickle.dump(self.meta_c_dict, mw)
        with open(meta_c, 'rb') as mr:
            self.meta_c_dict = pickle.load(mr)
        # hmdb 有kegg C号的化合物
        # meta_hmdb = 'Y:\\蛋白代谢事业部\\蛋白与代谢生信分析部\\06实验室研发项目\\sx_do_in_win\\代谢数据库\\meta_hmdb.pk'
        meta_hmdb = meta_hmdb_pk
        if not os.path.exists(meta_hmdb):
            # with open('Y:\\蛋白代谢事业部\\蛋白与代谢生信分析部\\06实验室研发项目\\sx_do_in_win\\代谢数据库\\hmdb_metabolites.detail.xls', 'rb') as hr:
            with open(meta_hmdb_pk_, 'rb') as hr:
                hmdb_df = pd.read_csv(hr, sep='\t', usecols=['accession', 'kegg_id'], index_col=0).dropna(how='any')
            self.meta_hmdb_dict = {ind: 1 for ind in hmdb_df.index}
            with open(meta_hmdb, 'wb') as mw:
                pickle.dump(self.meta_hmdb_dict, mw)
        with open(meta_hmdb, 'rb') as mr:
            self.meta_hmdb_dict = pickle.load(mr)
        # 读入外源性代谢物，这个量比较少，我就直接读入了
        # with open('Y:\\蛋白代谢事业部\\蛋白与代谢生信分析部\\06实验室研发项目\\sx_do_in_win\\代谢数据库\\extron_meta.txt') as em:
        with open(extron_meta) as em:
            self.extro_dict = {line.strip().split('\t')[-1]: 1 for line in em if line.strip()}

    #     进行分组字典建立与改名字典建立
        self.group2samples = dict()
        if os.path.exists(group_table):
            with open(group_table) as gr:
                df_group = pd.read_csv(gr, dtype={0: str, 1: str}, sep='\t')
            for group, sample in zip(df_group['group'], df_group['#sample']):
                if group not in self.group2samples:
                    self.group2samples[group] = list()
                if not sample in self.group2samples[group]:
                    self.group2samples[group].append(sample)
        else:
            print('警告：没有发现group文件')
        self.changename = dict()
        if os.path.exists(chname):
            with open(chname) as cr:
                chname_df = pd.read_csv(cr, dtype={0: str, 1: str}, sep='\t', header=None)
            self.changename = dict(zip(chname_df[0], chname_df[1]))
            # tmp = copy.copy(self.changename)
            # for c in tmp:
            #     if c.endswith('-n'):
            #         self.changename[c.replace('-n', '-p')] = self.changename[c]
            #     if c.endswith('-p'):
            #         self.changename[c.replace('-p', '-n')] = self.changename[c]
        else:
            print('警告：没有发现changename文件')
        if '%' in cutoff_percent:
            self.cutoff_percent = float(cutoff_percent.strip('%'))/100
        else:
            self.cutoff_percent = float(cutoff_percent)
        # 80%原则是指要有80%的非零值才保留，所以这边需要用1减一下
        self.cutoff_percent = 1 - self.cutoff_percent

        # 为了整合显红之前写的脚本，需要读入注释信息等，这边直接复制了显红的脚本，并稍微进行了修改
        # anno_info = 'Y:\\蛋白代谢事业部\\蛋白与代谢生信分析部\\06实验室研发项目\\sx_do_in_win\\代谢数据库\\anno_info.pk'
        anno_info =anno_info_pk #需要更改20220124
        # if not os.path.exists(anno_info):
        #     self.deal_anno_info(anno_info)
        # with open(anno_info, 'rb') as an:
        #     self.dic_name2cid, self.dic_name2cpd, self.dic_name2cas, self.dic_name2formula = pickle.load(an)
        #2022-01-24根据封一统需求更改注释脚本
        # anno_info = 'Y:\\蛋白代谢事业部\\蛋白与代谢生信分析部\\06实验室研发项目\\sx_do_in_win\\代谢数据库\\anno_info_new.pk'  #直接采用了显红新整理的配置注释的脚本
        if not os.path.exists(anno_info):
            if anno_info.endswith('_new.pk'):
                print('采用显红在2022年新整理的注释脚本')
                from get_pubid_info_dict import get_public_anno_info
                get_public_anno_info(anno_info,self.cid2detail,self.hmdb,self.lipids,smpdb_file)
            else:
                self.deal_anno_info(anno_info)
        with open(anno_info, 'rb') as an:
            self.dic_name2cid, self.dic_name2cpd, self.dic_name2cas, self.dic_name2formula = pickle.load(an)

        # 有老师特别关注的物质，老师不样删除
        self.free_f_d = dict()
        if os.path.exists(filter_free):
            with open(filter_free) as ff:
                for line in ff:
                    if not line.strip():
                        continue
                    line = line.strip().split('\t')[0]
                    self.free_f_d[line] = 1
                    self.free_f_d[self.get_stardard(line)] =1
        # print(self.free_f_d)
        self.filter_free = filter_free
        # 后来LC的项目都增加了内参，这个也需要针对内参问题处理一下
        self.inter_c_d = dict()
        if inter_compound and os.path.exists(inter_compound):
            with open(inter_compound) as ir:
                self.inter_c_d = {line.strip().split('\t')[0]: 1 for line in ir if line.strip()}
            self.free_f_d.update(self.inter_c_d)
            self.free_f_d['Internal standard'] = 1
            if not self.filter_free:
                self.filter_free = inter_compound
            else:
                with open('combined_free.list', 'w') as cw:
                    cw.write('\n'.join(self.free_f_d.keys()))
                self.filter_free = 'combined_free.list'
        # 美吉自建库之后，将ID对应到代谢物名称
        # MJ_db_p = pd.read_csv(r'Y:\蛋白代谢事业部\蛋白与代谢生信分析部\06实验室研发项目\sx_do_in_win\代谢数据库\MJDB_POS.20210621.v6.msp.xls', sep='\t', usecols=['ID', 'NAME']).drop_duplicates()
        # MJ_db_n = pd.read_csv(r'Y:\蛋白代谢事业部\蛋白与代谢生信分析部\06实验室研发项目\sx_do_in_win\代谢数据库\MJDB_NEG.20210621.v6.msp.xls', sep='\t', usecols=['ID', 'NAME']).drop_duplicates()
        # MJ_db_df = pd.concat([MJ_db_p, MJ_db_n])
        # MJ_db_df = pd.read_csv(r'Y:\蛋白代谢事业部\蛋白与代谢生信分析部\06实验室研发项目\sx_do_in_win\代谢数据库\MJDB_selfbuild.basic_info_v10.xls',sep='\t', usecols=['ID', 'NAME']).drop_duplicates()
        MJ_db_df = pd.read_csv(r'{}'.format(basic_info_v10),sep='\t', usecols=['ID', 'NAME']).drop_duplicates()
        self.MJ_db = dict(zip(MJ_db_df['ID'], MJ_db_df['NAME']))
        self.des2dbid = dict()
        os.chdir(os.path.dirname(self.pos_q))
        # 2020-10-09新增屏幕输出同时保存到某个文件的功能
        # try:
        #     import time
        #     from pplogger import Logger
        #     file = time.strftime("%Y-%m-%d", time.localtime()) + '.log.txt'
        #     sys.stdout = Logger(file, stream=sys.stdout)
        # except Exception as e:
        #     print(e)
        #     print('警告： 屏幕输出信息无法同时输出到文件')
        #     pass

    # 用来处理代谢物的名称
    def get_stardard(self, name):
        try:
            for alpha in self.dict_al2al:
                if alpha in name:
                    name = name.replace(alpha, self.dict_al2al[alpha])
        except:
            pass
        if ' / ' in name:
            name = name.split(' / ')[0]
        name_len = len(name)
        if name_len > 0:
            name_head = name[0]
            name_tail = name[1:name_len + 1]
            name = ''.join([name_head.upper(), name_tail])
        return name

    def deal_description(self, row):
        comp = str(row['Compound'])
        if comp in self.inter_c_d:
            return 'Internal standard'
        # des = row['Description']
        des = row['Description'].replace('"', "'") # 去掉双引号
        if des.startswith('MJ'):
            dbid = des
            des = self.MJ_db.get(des, des)
            self.des2dbid[des] = dbid
            self.des2dbid[self.get_stardard(des)] = dbid
        des = self.get_stardard(des)
        # 左旋氨基酸都改成右旋，没有旋的加L-
        if des.startswith('D-'):
            tmp = self.amino_acid.get(des.split('D-')[1])
            if tmp:
                des = 'L-' + des.split('D-')[1]
        if des.startswith('DL-'):
            tmp = self.amino_acid.get(des.split('DL-')[1])
            if tmp:
                des = 'L-' + des.split('DL-')[1]
        if des in self.amino_acid:
            des = 'L-' + des
        return des

    # 通过score筛选代谢物的函数，Fragmentation Score或Theoretical Fragmentation Score有一项达标就可以了
    def filter_by_scores(self, row, cutoff=30):
        des = row['Description']
        if self.free_f_d and (self.free_f_d.get(des) or self.free_f_d.get(des.lstrip('L-'))):
            return True
        score = row['Fragmentation Score']
        if score.__str__()[0].isdigit():
            if float(score) > cutoff:
                return True
        try:
            score = row['Theoretical Fragmentation Score']
        except:
            return False
        if score.__str__()[0].isdigit():
            if float(score) > cutoff:
                return True
        return False

    # 通过代谢物名称和ID等筛选代谢物，只有有C号的就留下，外源的不要，多肽不要，名字太长的不要
    def filter_by_description(self, row):
        des = row['Description']
        # 不让筛去的就不筛
        if self.free_f_d.get(des):
            return True
        # 判断是否是外源物，是的返回F
        if self.extro_dict.get(des):
            return False
        # 如果是HMDB的则判断是否有C号，有的马上返回T
        cid = row['Compound ID']
        if cid.__str__().startswith('HMDB'):
            if self.meta_hmdb_dict.get(cid):
                return True

        # 看description信息有没有C号，有的返回T
        if self.meta_c_dict.get(des):
            return True
        # 经历了上述判断，再看代谢物的长度，如果长度过长，大于80，则返回F
        if len(des) > self.max_len:
            return False
        # 如果是多肽则返回F，但是有些情况下客户会选择保留多肽
        if self.kegg_amino_acid == 'no':
            for i in des.split(' '):
                if len(i) != 3:
                    break
                if not self.peps.get(i):
                    break
            else:
                return False
        # 经历完所有判断的默认为T
        return True

    def too_small_is_zero(self, x, cut=4):
        if round(x, cut) > 0:
            return x
        return 0

    # 删除缺失值过多的代谢物
    def filter_toomuch_na_metas(self, row):
        comp = str(row['Metabolite'])
        if comp in self.free_f_d:
            print('%s没有参加na过滤'%comp)
            return True
        def count_zero(quants):
            count = 0
            for i in quants:
                if not float(i) > 0:
                    count += 1
            return count
        if not self.group2samples:
            quants = list()
            for n, ind in enumerate(row.index):
                if n > 12 and 'QC' not in ind:
                    quants.append(row[ind])
            if not len(quants):
                return True
            if float(count_zero(quants))/len(quants) < self.cutoff_percent:
                return True
            else:
                return False
        for group, samples in self.group2samples.items():
            if group == 'QC':
                continue
            quants = row.loc[samples].tolist()
            if float(count_zero(quants))/len(quants) < self.cutoff_percent:
                return True
        return False

    # 修改样品名字，返回新的columns列表
    def get_new_columns(self, column):
        if type(column) != list:
            column = column.tolist()
        for n, c in enumerate(column):
            tmp = self.changename.get(c)
            if tmp:
                column[n] = tmp
                continue
            if 'neg_' in c or 'pos_' in c or 'NEG_' in c or 'POS_' in c:
                tmp = '_'.join(c.split('_')[2:])
                column[n] = tmp
                tmp = self.changename.get(tmp)
                if tmp:
                    column[n] = tmp
        return column

    # 数据框根据description等信息去重，返回的是保留列的index
    def get_unique_index(self, df):
        unique_dict = dict()
        prefer_list = list()
        des_ = 'Description'
        df['Description'] = df['Description'].astype(str)
        if 'Theoretical Fragmentation Score' not in df.columns:
            des_ = 'Compound'
        for ind in df.index:
            des = str(df.loc[ind, des_])
            # 我也不知道为啥会有‘\n’--后来发现是俩数据框合并的时候没有选择ignore_index参数，失策失策
            des = des.split('\n')[0]
            if '   ' in des:
                des = des.split('  ')[1].strip()
            # 要保留的代谢需要跳过去重
            if des in self.free_f_d:
                prefer_list.append(ind)
                print(des+'跳过了去重')
                continue
            # 貌似不能区分大小写啊
            des = des.lower()
            score = df.loc[ind, 'Fragmentation Score'].__str__()
            if not score[0].isdigit():
                score = 0
            try:
                tscore = df.loc[ind, 'Theoretical Fragmentation Score'].__str__()
            except:
                tscore = '0'
            if not tscore[0].isdigit():
                tscore = 0
            ind_ = unique_dict.get(des)
            if not ind_:
                unique_dict[des] = ind
            else:
                score_ = df.loc[ind_, 'Fragmentation Score'].__str__()
                if not score_[0].isdigit():
                    score_ = 0
                try:
                    tscore_ = df.loc[ind_, 'Theoretical Fragmentation Score'].__str__()
                except:
                    tscore_ = '0'
                if not tscore_[0].isdigit():
                    tscore_ = 0
                if float(score_) < float(score):
                    unique_dict[des] = ind
                    continue
                if float(tscore_) < float(tscore):
                    unique_dict[des] = ind
        unique_dict = {v: k for k, v in unique_dict.items()}
        unique_index = list(unique_dict.keys())
        unique_index.sort()
        for ind in prefer_list:
            if ind not in unique_index:
                unique_index.insert(0, ind)
            else:
                unique_index.remome(ind)
                unique_index.insert(0, ind)
        return unique_index


    def get_unique_index_mini(self, df):
        unique_dict = dict()
        des_ = 'Description'
        df['Description'] = df['Description'].astype(str)
        for ind in df.index:
            des = str(df.loc[ind, des_])
            # 我也不知道为啥会有‘\n’--后来发现是俩数据框合并的时候没有选择ignore_index参数，失策失策
            des = des.split('\n')[0]
            if '   ' in des:
                des = des.split('  ')[1].strip()
            # 要保留的代谢需要跳过去重
            # 貌似不能区分大小写啊
            des = des.lower()
            score = df.loc[ind, 'Fragmentation Score'].__str__()
            if not score[0].isdigit():
                score = 0
            try:
                tscore = df.loc[ind, 'Theoretical Fragmentation Score'].__str__()
            except:
                tscore = '0'
            if not tscore[0].isdigit():
                tscore = 0
            ind_ = unique_dict.get(des)
            if not ind_:
                unique_dict[des] = ind
            else:
                score_ = df.loc[ind_, 'Fragmentation Score'].__str__()
                if not score_[0].isdigit():
                    score_ = 0
                try:
                    tscore_ = df.loc[ind_, 'Theoretical Fragmentation Score'].__str__()
                except:
                    tscore_ = '0'
                if not tscore_[0].isdigit():
                    tscore_ = 0
                if float(score_) < float(score):
                    unique_dict[des] = ind
                    continue
                if float(tscore_) < float(tscore):
                    unique_dict[des] = ind
        unique_dict = {v: k for k, v in unique_dict.items()}
        unique_index = list(unique_dict.keys())
        unique_index.sort()
        return unique_index


    # 处理feasureid
    def deal_feasureid(self, row, type_='pos'):
        try:
            feasure = row['Feature ID'].__str__()
        except:
            feasure = row.name.__str__()
        return type_ + '_' + feasure
        # return type_ + feasure

    # 给表格排序
    def get_index(self, ele, select_useful=False, columns=None):
        sorted_header =[
            "ID", "Metabolite", "m/z", "Retention time", "Mode",
            "Adducts", "Formula", "Fragmentation Score",
            "Theoretical Fragmentation Score", "Mass Error (ppm)",
            "KEGG Compound ID", "Library ID", "CAS ID"
        ]
        count_header = sum([v for k, v in self.group2samples.items() if k != 'QC'], [])
        sorted_header += count_header
        # 有的项目会有多余的没见过的列,所以给这个函数增加了一个功能
        if select_useful:
            useful_cols = list()
            for col in columns:
                if col in sorted_header or 'QC' in col:
                    useful_cols.append(col)
            return useful_cols
        try:
            return sorted_header.index(ele)
        except:
            return len(sorted_header) + 1

    # 运行逻辑控制步骤
    def run(self):
        # 现在需要先判断这是不是未经过生产处理的原始文件
        with open(self.pos_q) as pr:
            headers = pr.readline().strip()
            if headers.startswith(','):
                from predeal_LC_raw import PredealLcRaw
                raw_predeal = PredealLcRaw(self.pos_i, self.pos_q, self.neg_i, self.neg_q, self.output, filter_free=self.filter_free,keep_hmdb=self.keep_hmdb, p_object=self)
                self.pos_i, self.pos_q, self.neg_i, self.neg_q = raw_predeal.run()

        # 读入各个文件并改样本名
        print('开始读入并处理各个文件')
        # 确定一下文件的编码
        def get_encoding(file):
            # 二进制方式读取，获取字节数据，检测类型
            with open(file, 'rb') as f:
                return chardet.detect(f.read())['encoding']
        encoding = 'GB18030'
        if get_encoding(self.pos_i) == 'utf-8':
            encoding = 'utf-8'

        drop_c = ['Neutral mass (Da)', 'Charge', 'Isotope Distribution']
        # if encoding != 'utf-8':
        #     pos_i_df = pd.read_csv(self.pos_i, sep=',', quoting=0, quotechar='"', encoding=encoding)
        #     pos_i_df.to_csv('tmp.xls', sep='\t', quoting=0, quotechar='"', index=False, header=True,
        #                     encoding='utf-8')
        #     pos_i_df = pd.read_csv('tmp.xls', sep='\t', quoting=3, encoding='utf-8')
        # else:
        #     pos_i_df = pd.read_csv(self.pos_i, sep=',', quoting=0, quotechar='"', encoding=encoding).fillna(0).replace('_', '-')
        with open(self.pos_i, 'rb') as pir:
            pos_i_df = pd.read_csv(pir, sep=',', quoting=0, quotechar='"', encoding=encoding).fillna(0).replace('_', '-')
        pos_i_df.columns = self.get_new_columns(pos_i_df.columns)
        with open(self.pos_q, 'rb') as pqr:
            pos_q_df = pd.read_csv(pqr, sep=',', quoting=3).fillna(0)
        # 偶尔给的数据会有多个空列，产生Unname
        pos_q_df = pos_q_df.loc[:, ~pos_q_df.columns.str.contains('^Unnamed')]
        pos_q_df = pos_q_df.drop(columns=drop_c)
        pos_q_df['Feature ID'] = pos_q_df.apply(self.deal_feasureid, type_='pos', axis=1)
        pos_q_df.columns = self.get_new_columns(pos_q_df.columns)
        #pos_q_df.to_csv('pos_q_df_raw.xls', sep='\t', quoting=3, index=False, header=True)
        # if encoding != 'utf-8':
        #     neg_i_df = pd.read_csv(self.neg_i, sep=',', quoting=0, quotechar='"', encoding=encoding)
        #     neg_i_df.to_csv('tmp.xls', sep='\t', quoting=0, quotechar='"', index=False, header=True,
        #                     encoding='utf-8')
        #     neg_i_df = pd.read_csv('tmp.xls', sep='\t', quoting=3, encoding='utf-8')
        # else:
        #     neg_i_df = pd.read_csv(self.neg_i, sep=',', quoting=0, quotechar='"', encoding=encoding).fillna(0).replace('_', '-')
        with open(self.neg_i, 'rb') as nir:
            neg_i_df = pd.read_csv(nir, sep=',', quoting=0, quotechar='"', encoding=encoding).fillna(0).replace('_', '-')
        neg_i_df.columns = self.get_new_columns(neg_i_df.columns)
        with open(self.neg_q, 'rb') as nqr:
            neg_q_df = pd.read_csv(nqr, sep=',', quoting=3).fillna(0)
        neg_q_df = neg_q_df.loc[:, ~neg_q_df.columns.str.contains('^Unnamed')]
        neg_q_df = neg_q_df.drop(columns=drop_c)
        neg_q_df['Feature ID'] = neg_q_df.apply(self.deal_feasureid, type_='neg', axis=1)
        neg_q_df.columns = self.get_new_columns(neg_q_df.columns)


        #neg_q_df.to_csv('neg_q_df_raw.xls', sep='\t', quoting=3, index=False, header=True)
        # 定性表通过最小分值筛选

        # 生成人员有时候会忘了删除4列，就很惆怅，害我们多写一步
        may_drop_c = ['Chromatographic peak width (min)', 'Identifications', 'Maximum Abundance', 'Minimum CV%']
        for col in may_drop_c:
            if col in pos_q_df.columns:
                pos_q_df = pos_q_df.drop(columns=col)
            if col in neg_q_df.columns:
                neg_q_df = neg_q_df.drop(columns=col)

        min_score = self.min_score
        if self.sure_score:
            min_score = self.sure_score
        print('以最小得分%s筛选定性数据' % min_score)
        pos_i_df = pos_i_df[pos_i_df.apply(self.filter_by_scores, cutoff=min_score, axis=1)]
        print('pos以得分%s筛选，筛选完之后还有%s行' % (str(min_score), str(pos_i_df.shape[0])))
        neg_i_df = neg_i_df[neg_i_df.apply(self.filter_by_scores, cutoff=min_score, axis=1)]
        print('neg以得分%s筛选，筛选完之后还有%s行' % (str(min_score), str(neg_i_df.shape[0])))

        # 处理定性表
        # print('正在通过Description筛选定性数据')
        # pos_i_df['Description'] = pos_i_df.apply(self.deal_description, axis=1)
        # pos_i_df = pos_i_df[pos_i_df.apply(self.filter_by_description, axis=1)]
        # neg_i_df['Description'] = neg_i_df.apply(self.deal_description, axis=1)
        # neg_i_df = neg_i_df[neg_i_df.apply(self.filter_by_description, axis=1)]

        # 加入注释信息
        print('正在添加注释信息')
        anno_cols = ["KEGG Compound ID", "Library ID", "CAS ID"]
        tmp_s = pos_i_df.apply(self.get_anno_info, axis=1)
        for n, i in enumerate(anno_cols):
            pos_i_df[i] = tmp_s.apply(lambda x: x[n])
        tmp_s = neg_i_df.apply(self.get_anno_info, axis=1)
        for n, i in enumerate(anno_cols):
            neg_i_df[i] = tmp_s.apply(lambda x: x[n])
        #neg_i_df.to_csv('neg_i_df_raw.xls', sep='\t', quoting=3, index=False, header=True)
        # pos和neg合并后去重
        print('正在进行去重操作')
        pos_i_df['TYPE'] = 'POS'
        neg_i_df['TYPE'] = 'NEG'
        pos_neg = pd.concat([pos_i_df, neg_i_df], ignore_index=True)
        uniq_ = self.get_unique_index(pos_neg)
        pos_neg = pos_neg.loc[uniq_,]

        if self.free_f_d:
            print('检测到有关注物质，多去重一次，避免关注物质在正负离子表都有')
            uniq_ = self.get_unique_index_mini(pos_neg)
            pos_neg = pos_neg.loc[uniq_,]

        # 去重之后需要控制pos和neg加起来不超过700，但是筛选用的值不超过50
        # 这几句注释掉是换了筛选策略，要设置三个档，30 40 50 得分筛，所以这个逻辑要换一下了
        # pn_n = self.min_score
        # while pos_neg.shape[0] > self.max_row:
        #     pn_n += 0.1
        #     if pn_n > self.max_score:
        #         break
        #     pos_neg = pos_neg[pos_neg.apply(self.filter_by_scores, cutoff=pn_n, axis=1)]
        if not self.sure_score:
            if pos_neg.shape[0] > self.max_row:
                print('以得分30筛选，筛选完之后还有%s行' % str(pos_neg.shape[0]))
                pos_neg_40 = pos_neg[pos_neg.apply(self.filter_by_scores, cutoff=40, axis=1)]
                print('以得分40筛选，筛选完之后还有%s行'%str(pos_neg_40.shape[0]))
                if self.min_row <= pos_neg_40.shape[0] <= self.max_row:
                    pos_neg = pos_neg_40
                if pos_neg_40.shape[0] > self.max_row:
                    pos_neg = pos_neg_40
                    pos_neg_50 = pos_neg_40[pos_neg_40.apply(self.filter_by_scores, cutoff=50, axis=1)]
                    print('以得分50筛选，筛选完之后还有%s行' % str(pos_neg_50.shape[0]))
                    if pos_neg_50.shape[0] < self.min_row:
                        pass
                    else:
                        pos_neg = pos_neg_50
        else:
            pos_neg = pos_neg[pos_neg.apply(self.filter_by_scores, cutoff=self.sure_score, axis=1)]

        pos_i_df = pos_neg[pos_neg['TYPE'] == 'POS']
        neg_i_df = pos_neg[pos_neg['TYPE'] == 'NEG']


        with open(os.path.join(self.output, 'pos_identification_raw.xls'), 'wb') as prw:
            tmp_f = StringIO()
            pos_i_df.to_csv(tmp_f, sep='\t', quoting=3, index=False, header=True, encoding='utf-8')
            prw.write(tmp_f.getvalue().encode('utf8'))
            tmp_f.close()
        with open(os.path.join(self.output, 'neg_identification_raw.xls'), 'wb') as nrw:
            tmp_f = StringIO()
            neg_i_df.to_csv(tmp_f, sep='\t', quoting=3, index=False, header=True,
                        encoding='utf-8')
            nrw.write(tmp_f.getvalue().encode('utf8'))
            tmp_f.close()

        if not 'Theoretical Fragmentation Score' in pos_i_df.columns.tolist():
            pos_i_df['Theoretical Fragmentation Score'] = '_'
            neg_i_df['Theoretical Fragmentation Score'] = '_'
        # 讲定性信息整合到定量信息上
        select_cols = ['Compound', 'Adducts', 'Formula', 'Fragmentation Score',
                       'Theoretical Fragmentation Score', 'Mass Error (ppm)', 'Description', "KEGG Compound ID", "Library ID", "CAS ID"]

        pos_i_df = pos_i_df[select_cols]
        neg_i_df = neg_i_df[select_cols]
        rename_dict = {
            'Description': 'Metabolite',
            'Feature ID': 'ID',
            # 'Mass Error (ppm)': 'Mass Error',
            'Retention time (min)': 'Retention time'
        }

        print('正在合并定性定量数据表')
        pos_i_q = pd.merge(pos_q_df, pos_i_df, on='Compound', how='left').fillna('-')
        pos_i_q = pos_i_q.rename(columns=rename_dict)
        pos_i_q['Mode'] = 'pos'
        c_pos = pos_i_q['Compound']
        pos_i_q = pos_i_q.drop(columns=['Compound'])
        pos_i_q_cols = pos_i_q.columns.tolist()
        pos_i_q_cols.sort(key=self.get_index)
        pos_i_q = pos_i_q[pos_i_q_cols]
        if self.group2samples:
            pos_i_q = pos_i_q[self.get_index(None, select_useful=True, columns=pos_i_q.columns)]

        neg_i_q = pd.merge(neg_q_df, neg_i_df, on='Compound', how='left').fillna('-')
        neg_i_q = neg_i_q.rename(columns=rename_dict)
        neg_i_q['Mode'] = 'neg'
        c_neg = neg_i_q['Compound']
        neg_i_q = neg_i_q.drop(columns=['Compound'])
        neg_i_q_cols = neg_i_q.columns.tolist()
        neg_i_q_cols.sort(key=self.get_index)
        neg_i_q = neg_i_q[neg_i_q_cols]
        if self.group2samples:
            neg_i_q = neg_i_q[self.get_index(None, select_useful=True, columns=neg_i_q.columns)]

        if self.free_f_d:
            print('检测到有关注物质，劣币驱逐良币，将跟关注物种一个compund id的物质去除')
            pos_i_q_new = pd.DataFrame(columns=pos_i_q.columns)
            for com, df in pos_i_q.groupby('ID'):
                if df.shape[0] > 1:
                    tmp_df = df[df['Metabolite'].isin(self.free_f_d)]
                    # ids = tmp_df['ID'].tolist()
                    # if len(ids) > 1:
                    #     tmp_df['ID'] = [id_ + '_' + str(n + 1) for n, id_ in enumerate(ids)]
                    pos_i_q_new = pos_i_q_new.append(tmp_df.iloc[0, :])
                else:
                    pos_i_q_new = pos_i_q_new.append(df)
            pos_i_q = pos_i_q_new
            neg_i_q_new = pd.DataFrame(columns=neg_i_q.columns)
            for com, df in neg_i_q.groupby('ID'):
                if df.shape[0] > 1:
                    tmp_df = df[df['Metabolite'].isin(self.free_f_d)]
                    # ids = tmp_df['ID'].tolist()
                    # if len(ids) > 1:
                    #     tmp_df['ID'] = [id_ + '_' + str(n + 1) for n, id_ in enumerate(ids)]
                    neg_i_q_new = neg_i_q_new.append(tmp_df.iloc[0, :])
                else:
                    neg_i_q_new = neg_i_q_new.append(df)
            neg_i_q = neg_i_q_new

        with open(os.path.join(self.output, 'pos_measurement_raw.xls'), 'wb') as prw:
            tmp_f = StringIO()
            pos_i_q.to_csv(tmp_f, sep='\t', quoting=3, index=False, header=True)
            prw.write(tmp_f.getvalue().encode('utf8'))
            tmp_f.close()
        with open(os.path.join(self.output, 'neg_measurement_raw.xls'), 'wb') as nrw:
            tmp_f = StringIO()
            neg_i_q.to_csv(tmp_f, sep='\t', quoting=3, index=False, header=True)
            nrw.write(tmp_f.getvalue().encode('utf8'))
            tmp_f.close()


        # 将最后的表做一个筛选工作
        # 在筛选之前需要把特别特别小的值，强行变成0
        for n, col in enumerate(pos_i_q.columns):
            if n > 12 and 'QC' not in col:
                pos_i_q[col] = pos_i_q[col].apply(self.too_small_is_zero)
        for n, col in enumerate(neg_i_q.columns):
            if n > 12 and 'QC' not in col:
                neg_i_q[col] = neg_i_q[col].apply(self.too_small_is_zero)
        print('正在通过每组阙值进行筛选')
        pos_i_q = pos_i_q[pos_i_q.apply(self.filter_toomuch_na_metas, axis=1)]
        neg_i_q = neg_i_q[neg_i_q.apply(self.filter_toomuch_na_metas, axis=1)]
        with open(os.path.join(self.output, 'pos_measurement_nodelcid.xls'), 'wb') as pmw:
            tmp_f = StringIO()
            pos_i_q.to_csv(tmp_f, sep='\t', quoting=3, index=False, header=True, encoding='utf-8')
            pmw.write(tmp_f.getvalue().encode('utf8'))
            tmp_f.close()
        with open(os.path.join(self.output, 'neg_measurement_nodelcid.xls'), 'wb') as nmw:
            tmp_f = StringIO()
            neg_i_q.to_csv(tmp_f, sep='\t', quoting=3, index=False, header=True, encoding='utf-8')
            nmw.write(tmp_f.getvalue().encode('utf8'))
            tmp_f.close()
        pos_i_q.insert(0, 'Compound', c_pos)
        with open(os.path.join(self.output, 'pos_measurement_nodelcid.xlsx'), 'wb') as pmw:
            pos_i_q.to_excel(pmw, index=False, header=True, encoding='utf-8')
        neg_i_q.insert(0, 'Compound', c_neg)
        with open(os.path.join(self.output, 'neg_measurement_nodelcid.xlsx'), 'wb') as nmw:
            neg_i_q.to_excel(nmw, index=False, header=True, encoding='utf-8')
        print('开始处理pos结果表里cid有多个的情况')
        try:
            self.adjust_mCID(os.path.join(self.output, 'pos_measurement_nodelcid.xls'), os.path.join(self.output, 'pos_measurement.xls'))
        except Exception as e:
            print('adjust_mCID出错,%s', e)
            print('pos不需要处理')
            try:
                os.remove(os.path.join(self.output, 'pos_measurement.xls'))
                os.remove(os.path.join(self.output, 'pos_measurement.xlsx'))
            except Exception:
                pass
            os.rename(os.path.join(self.output, 'pos_measurement_nodelcid.xls'), os.path.join(self.output, 'pos_measurement.xls'))
            os.rename(os.path.join(self.output, 'pos_measurement_nodelcid.xlsx'), os.path.join(self.output, 'pos_measurement.xlsx'))
        print('开始处理neg结果表里cid有多个的情况')
        try:
            self.adjust_mCID(os.path.join(self.output, 'neg_measurement_nodelcid.xls'), os.path.join(self.output, 'neg_measurement.xls'))
        except Exception as e:
            print('adjust_mCID出错,%s', e)
            print('neg不需要处理')
            try:
                os.remove(os.path.join(self.output, 'neg_measurement.xls'))
                os.remove(os.path.join(self.output, 'neg_measurement.xlsx'))
            except Exception:
                pass
            os.rename(os.path.join(self.output, 'neg_measurement_nodelcid.xls'), os.path.join(self.output, 'neg_measurement.xls'))
            os.rename(os.path.join(self.output, 'neg_measurement_nodelcid.xlsx'), os.path.join(self.output, 'neg_measurement.xlsx'))

    def deal_anno_info(self, out_file):
        from collections import defaultdict
        dic_name2cid = defaultdict(list)
        dic_name2cpd = defaultdict(list)
        dic_name2cas = defaultdict(list)
        dic_name2formula = defaultdict(list)
        with open(self.kegg, 'r')as KO2name_r:
            for line in KO2name_r.readlines():
                if '\t' in line:
                    items = line.strip().split('\t')
                    CID = items[0].split(':')[1]
                    names = items[1].split('; ')
                    for name in names:
                        name = name.lower()
                        dic_name2cid[name].append(CID)
                        dic_name2cid[self.get_stardard(name)].append(CID)

        with open(self.hmdb, 'rb')as hmdb_r:
            for line in hmdb_r.readlines():
                line = line.decode('utf8')
                if '\t' in line:
                    items = line.strip('\n').split('\t')
                    CID = items[8]
                    CAS = items[7]
                    CPD = items[0]
                    FORMULA = items[2]
                    names = items[1].split(';') + items[5].split(';') + items[6].split(';')
                    names_ = items[1].split(';') + items[5].split(';') + items[-1].split(';') + items[6].split(';')
                    names_cid = list(filter(None, names_))
                    for name in names:
                        name = name.lower()
                        if CAS != "_":
                            dic_name2cas[name].append(CAS)
                            dic_name2cas[self.get_stardard(name)].append(CAS)
                        dic_name2formula[name].append(FORMULA)
                        dic_name2formula[self.get_stardard(name)].append(FORMULA)
                        # dic_name2cpd[name].append(CPD)
                        # dic_name2cpd[self.get_stardard(name)].append(CPD)
                    for name_cid in names_cid:
                        name_cid = name_cid.lower()
                        if CID != "_":
                            dic_name2cid[name_cid].append(CID)
                            dic_name2cid[self.get_stardard(name_cid)].append(CID)
                        dic_name2cpd[name_cid].append(CPD)
                        dic_name2cpd[self.get_stardard(name_cid)].append(CPD)

        with open(self.lipids, 'rb')as lipids_r:
            for line in lipids_r.readlines():
                line = line.decode('utf8')
                if '\t' in line:
                    items = line.strip('\n').split('\t')
                    CPD = items[0]
                    CID = items[3]
                    FORMULA = items[8]
                    names = items[1].split(';') + items[2].split(';') + items[-1].split('; ')
                    for name in names:
                        if name in ['-', '_']:
                            continue
                        name = name.lower()
                        if CID != "_":
                            dic_name2cid[name].append(CID)
                            dic_name2cid[self.get_stardard(name)].append(CID)
                        dic_name2cpd[name].append(CPD)
                        dic_name2cpd[self.get_stardard(name)].append(CPD)
                        dic_name2formula[name].append(FORMULA)
                        dic_name2formula[self.get_stardard(name)].append(FORMULA)
        with open(out_file, 'wb') as ow:
            pickle.dump([dic_name2cid, dic_name2cpd, dic_name2cas, dic_name2formula], ow)
    # 忽略大小写匹配--戴显红修改
    def get_anno_info(self, row):
        def get_value(des, dict_):
            des = des.lower()
            tmp = dict_.get(des)
            if tmp:
                return ';'.join(list(set(tmp)))
            des = des.replace(u'(\xb1)', u'(\xc2\xb1)')
            des = u'%s'%des
            tmp = dict_.get(des)
            if tmp:
                return ';'.join(list(set(tmp)))
            return '-'

        des = row['Description']
        cid = get_value(des, self.dic_name2cid)
        # 氨基酸全部为正弦的cid
        if des.startswith('L-'):
            tmp = self.amino_acid.get(des.split('L-')[1])
            if tmp:
                cid = tmp
                # print(des+':'+cid)
        cpd = get_value(des, self.dic_name2cpd)
        if des in self.des2dbid:
            if cpd == '-':
                cpd = ''
            cpd = cpd.split(';')
            cpd.append(self.des2dbid[des])
            cpd = ';'.join(cpd)
        cas = get_value(des, self.dic_name2cas)
        return cid, cpd, cas


    def adjust_mCID(self, file, out):
        # cid2name = r'Y:\蛋白代谢事业部\蛋白与代谢生信分析部\06实验室研发项目\sx_do_in_win\代谢数据库\CID2NAME.20210409.v98.0.txt'
        # cid2map = r'Y:\蛋白代谢事业部\蛋白与代谢生信分析部\06实验室研发项目\sx_do_in_win\代谢数据库\Compound_Pathway.list'
        # cid2detail = r'Y:\蛋白代谢事业部\蛋白与代谢生信分析部\06实验室研发项目\sx_do_in_win\代谢数据库\KEGG_compound.20210526_addHMDB.xls'
        cid2name = self.cid2name
        cid2map = self.cid2map
        cid2detail = self.cid2detail
        need_deal = DealMultiCid(file, cid2name, cid2map, cid2detail, out)
        need_deal.deal_each_meta()

# def __init__(self, pos_i, pos_q, neg_i, neg_q, max_row, min_score, max_len, group_table, chname, cutoff_percent)
# kegg_amino_acid = 'no'
# kegg = 'Y:\\蛋白代谢事业部\\蛋白与代谢生信分析部\\06实验室研发项目\\sx_do_in_win\\代谢数据库\\CID2name.txt'
# kegg = r'Y:\蛋白代谢事业部\蛋白与代谢生信分析部\06实验室研发项目\sx_do_in_win\代谢数据库\CID2NAME.20210409.v98.0.txt'
# hmdb = 'Y:\\蛋白代谢事业部\\蛋白与代谢生信分析部\\06实验室研发项目\\sx_do_in_win\\代谢数据库\\hmdb_metabolites.detail.xls'
# lipids = 'Y:\\蛋白代谢事业部\\蛋白与代谢生信分析部\\06实验室研发项目\\sx_do_in_win\\代谢数据库\\LIPIDMAPS.detail.xls'
# epa = 'Y:\\蛋白代谢事业部\\蛋白与代谢生信分析部\\06实验室研发项目\\sx_do_in_win\\代谢数据库\\EPA.detail.xls'
# parser = argparse.ArgumentParser(description="处理LC代谢的下机数据脚本")
#
# parser.add_argument("-kegg_amino_acid", choices=['yes', 'no'], type=str, default=kegg_amino_acid, help='有些客户会选择保留多肽')
# parser.add_argument("-kegg",
#                     default='Y:\\蛋白代谢事业部\\蛋白与代谢生信分析部\\06实验室研发项目\\sx_do_in_win\\代谢数据库\\CID2name.txt',
#                     help="input file as KO2name.txt")
# parser.add_argument("-hmdb",
#                     default='Y:\\蛋白代谢事业部\\蛋白与代谢生信分析部\\06实验室研发项目\\sx_do_in_win\\代谢数据库\\hmdb_metabolites.detail.xls',
#                     help="input file as class.hmdb_metabolites.xls")
# parser.add_argument("-lipids",
#                     default='Y:\\蛋白代谢事业部\\蛋白与代谢生信分析部\\06实验室研发项目\\sx_do_in_win\\代谢数据库\\LIPIDMAPS.detail.xls',
#                     help="input file as LIPIDMAPS.detail.xls")
# parser.add_argument("-epa",
#                     default='Y:\\蛋白代谢事业部\\蛋白与代谢生信分析部\\06实验室研发项目\\sx_do_in_win\\代谢数据库\\EPA.detail.xls',
#                     help="input file as EPA.detail.xls")
#
# args = parser.parse_args()


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description="")
    parser.add_argument("-pos_i", "--pos_i", dest="pos_i", required=True, type=str, help="")
    parser.add_argument("-pos_q", "--pos_q", dest="pos_q", required=True, type=str, help="")
    parser.add_argument("-neg_i", "--neg_i",dest="neg_i", required=True, type=str, help="")
    parser.add_argument("-neg_q", "--neg_q",dest="neg_q", required=True, type=str, help="")
    parser.add_argument("-min_row", "--min_row",dest="min_row", required=True, type=int, help="")
    parser.add_argument("-max_row", "--max_row",dest="max_row", required=True, type=int, help="")
    parser.add_argument("-min_score", "--min_score",dest="min_score", required=True, type=int, help="")
    parser.add_argument("-max_score", "--max_score",dest="max_score", required=True, type=int, help="")
    parser.add_argument("-sure_score", "--sure_score",dest="sure_score", required=True, type=int, help="")
    parser.add_argument("-max_len", "--max_len",dest="max_len", required=True, type=int, help="")
    parser.add_argument("-group_table", "--group_table",dest="group_table", required=True, type=str, help="")
    parser.add_argument("-chname", "--chname",dest="chname", required=True, type=str, help="")
    parser.add_argument("-cutoff_percent", "--cutoff_percent",dest="cutoff_percent", required=True, type=str, help="")
    parser.add_argument("-filter_free", "--filter_free",dest="filter_free", required=True, type=str, help="")
    parser.add_argument("-database_dir", "--database_dir",dest="database_dir", required=True, type=str, help="")
    parser.add_argument("-output", "--output",dest="output", required=True, type=str, help="")
    parser.add_argument("-kegg_amino_acid", "--kegg_amino_acid",dest="kegg_amino_acid", required=True, type=str, help="")
    args = parser.parse_args()
    predeal = PredealLcMatrix(args.pos_i, args.pos_q, args.neg_i, args.neg_q,\
        args.min_row, args.max_row, args.min_score, args.max_score, args.sure_score,\
        args.max_len, args.group_table, args.chname, args.cutoff_percent, args.filter_free,\
        None, args.database_dir, args.output,kegg_amino_acid=args.kegg_amino_acid)
    predeal.run()
