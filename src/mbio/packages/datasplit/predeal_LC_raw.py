#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""
@time    : 2019/9/25 14:20
@file    : predeal_LC_raw.py
@author  : yitong.feng
@contact: yitong.feng@majorbio.com
"""


import argparse
import os
import pandas as pd
import numpy as np
import re
from flower_print import print_flower_letter
from io import StringIO




class PredealLcRaw(object):
    def __init__(self, pos_i, pos_q, neg_i, neg_q, output_dir, use_col=0, filter_free=None,keep_hmdb='yes', p_object=None):
        # 传入四个表格和最大行数与最小得分
        self.pos_i = os.path.abspath(pos_i)
        self.pos_q = os.path.abspath(pos_q)
        self.neg_i = os.path.abspath(neg_i)
        self.neg_q = os.path.abspath(neg_q)
        self.output_dir = output_dir

        self.use_col_pos = use_col
        if not use_col:
            self.use_col_pos = self.get_use_col(self.pos_q)
        self.use_col_neg = use_col
        if not use_col:
            self.use_col_neg = self.get_use_col(self.neg_q)
        self.filter_free = dict()
        if filter_free and os.path.exists(filter_free):
            with open(filter_free) as fr:
                self.filter_free = {line.strip().split('\t')[0]: 1 for line in fr if line.strip()}
        self.keep_hmdb = keep_hmdb
        self.p_object = p_object

    # 得到定量表最后一列的index
    def get_use_col(self, table):
        with open(table) as pq:
            header = pq.readline().strip('\n').split(',')
            for n, h in enumerate(header):
                if 'Raw abundance' in h:
                    return n
            print('定量文件好像是有点问题，没有Raw abundance这一列')
            return n

    # 数据框根据Compound等信息去重，返回的是保留列的index
    def get_unique_index(self, df, use='Compound'):
        unique_dict = dict()
        perfer_dict = dict()
        des_ = use
        df['Description'] = df['Description'].astype(str)
        perfer_list = list()
        perfered_des = dict()
        for ind in df.index:
            des = str(df.loc[ind, des_])
            if df.loc[ind, 'Description'].lower() in perfered_des:
                continue
            com_id = str(df.loc[ind, 'Compound ID'])
            # 还得忽略一下大小写
            des = des.lower()
            score = df.loc[ind, 'Fragmentation Score']
            try:
                score = float(score)
            except:
                score = 0
            # 把判断是否是感兴趣代谢物放在了这边，减少计算，节省时间
            if self.filter_free and df.loc[ind, 'Compound'] in self.filter_free and df.loc[ind, 'Description'] in self.filter_free:
                print(df.loc[ind, 'Description']+'会保留')
                score = 103
                perfer = 2
                unique_dict[des] = ind
                perfer_dict[des] = perfer
                if not des in perfered_des:
                    perfer_list.append(ind)
                    perfered_des[df.loc[ind, 'Description'].lower()] = 1
                continue
            if com_id.startswith('HMDB') or score < 30:
                perfer = 0
            else:
                perfer = 1
            # 有时候会根据desciption留下一部分
            if self.filter_free and df.loc[ind, 'Description'] in self.filter_free:
                # print(df.loc[ind, 'Description']+'会保留')
                score = 101
                perfer = 1
            # 用compound会更精确，所以设置得分高一点
            if self.filter_free and df.loc[ind, 'Compound'] in self.filter_free:
                # print(df.loc[ind, 'Description']+'会保留')
                score = 102
                perfer = 1
            if not des in perfer_dict:
                perfer_dict[des] = perfer
            # ind_ = unique_dict.get(des)
            ind_ = None
            try:
                ind_ = unique_dict[des]
            except:
                pass
            if not ind_:
                unique_dict[des] = ind
            else:
                try:
                    perfer_ = perfer_dict[des]
                except:
                    perfer_ = 0
                if perfer_ == 2:
                    continue
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
        unique_index = list(set(unique_index + perfer_list))
        return unique_index

    # 将Fragmentation Score根据是不是HMDB分成两个
    def deal_score(self, row):
        com_id = row['Compound ID'].__str__()
        score = row['Fragmentation Score']
        try:
            score = float(score)
        except:
            score = 0
        if com_id.startswith('HMDB'):
            return -np.inf, score
        return score, -np.inf

    # 根据QC的值，将所有QC值都不足10的去掉
    def filter_qc(self, row, qc_cols, cutoff=10):
        fea = row['Compound']
        if self.filter_free.get(fea):
            return True
        if not qc_cols:
            return True
        qcs = row.loc[qc_cols]
        for qc in qcs:
            try:
                qc = float(qc)
            except:
                qc = 0
            if qc > cutoff:
                return True
        # print('%s filtered by QC'%row['Compound'])
        return False

    # 判断定量数据框里是rsd小于30的代谢物含量是否在70以上
    def judge_rsd(self, df, qc_cols, cutoff=0.3):
        def cal_rsd(row):
            qcs = row.loc[qc_cols]
            # rsd = round(qcs.std(), 5)/round(qcs.mean(), 5)
            rsd = qcs.std()/qcs.mean()
            return rsd
        rsd_ = df.apply(cal_rsd, axis=1)
        stard = rsd_[rsd_ < cutoff]
        if not rsd_.count():
            return 0
        return float(stard.count())/rsd_.count()


    # 运行逻辑控制步骤
    def run(self):
        # 读入各个文件并改样本名
        print('本项目判断为是从生产部门直接得到的数据\n现在直接从最原始的文件开始处理')
        # 处理定量文件，添加一列Feature ID,然后还要删除四列
        drop_c = ['Chromatographic peak width (min)', 'Identifications', 'Maximum Abundance', 'Minimum CV%']
        with open(self.pos_q, 'rb') as pqr:
            pos_q_df = pd.read_csv(pqr, sep=',', quoting=0, quotechar='"',
                                usecols=range(self.use_col_pos), skiprows=[0, 1])
        # pos_q_df = pos_q_df.drop(columns=drop_c)
        qcbfs = list()
        for c in pos_q_df.columns.tolist():
            if 'QCBF' in c:
                qcbfs.append(c)
        pos_q_df = pos_q_df.drop(columns=drop_c + qcbfs)
        # 实验室同事有时候还能给多搞出来五列，真的是没办法
        may_col = ["Anova (p)", "q Value", "Max Fold Change", "Highest Mean", "Lowest Mean"]
        for mc in may_col:
            if mc in pos_q_df.columns.tolist():
                pos_q_df = pos_q_df.drop(columns=[mc])
        re_qc = re.compile('QC\d')
        qc_cols = [col for col in pos_q_df.columns if re_qc.search(col)]
        pos_q_df = pos_q_df[pos_q_df.apply(self.filter_qc, qc_cols=qc_cols, axis=1)]
        rsd_v = self.judge_rsd(pos_q_df, qc_cols)
        if rsd_v < 0.7:
            print('pos定量表qc的rsd小于30的不足0.7，qc用20筛选')
            pos_q_df = pos_q_df[pos_q_df.apply(self.filter_qc, qc_cols=qc_cols, cutoff=20, axis=1)]
            rsd_v = self.judge_rsd(pos_q_df, qc_cols)
            if rsd_v < 0.7:
                if pos_q_df.shape[0] > 20000:
                    print('pos的总峰值超过20000，尝试qc用25筛选')
                    pos_q_df = pos_q_df[pos_q_df.apply(self.filter_qc, qc_cols=qc_cols, cutoff=25, axis=1)]
                    rsd_v = self.judge_rsd(pos_q_df, qc_cols)
                    if rsd_v < 0.7:
                        print_flower_letter('warnning', 'white', 'red', 'anti')
                        print('警告：pos定量表qc用25筛选后的rsd依旧不足0.7  现为%s' % str(rsd_v))
                else:
                    print_flower_letter('warnning', 'white', 'red', 'anti')
                    print('警告：pos定量表qc用20筛选后的rsd依旧不足0.7  现为%s' % str(rsd_v))
            if rsd_v < 0.7:
                # 尝试看看是否删除某列qc后rsd会达标
                for qc in qc_cols:
                    tmp = qc_cols[::]
                    tmp.remove(qc)
                    rsd_tmp = self.judge_rsd(pos_q_df, tmp)
                    if rsd_tmp >= 0.7:
                        print('你可以尝试删除%s' % qc)
                    print('去掉%s后，pos表的rsd为%s' % (qc, rsd_tmp))
        # pos_q_df.insert(0, 'Feature ID', pos_q_df.index + 1)
        print('本项目pos定量表的rsd为%s' % str(rsd_v))
        pos_q_df.insert(0, 'Feature ID', np.arange(0, pos_q_df.shape[0]) + 1)
        with open(self.neg_q, 'rb') as nqr:
            neg_q_df = pd.read_csv(nqr, sep=',', quoting=0, quotechar='"',
                               usecols=range(self.use_col_neg), skiprows=[0, 1])
        # neg_q_df = neg_q_df.drop(columns=drop_c)
        qcbfs = list()
        for c in neg_q_df.columns.tolist():
            if 'QCBF' in c:
                qcbfs.append(c)
        neg_q_df = neg_q_df.drop(columns=drop_c + qcbfs)
        for mc in may_col:
            if mc in neg_q_df.columns.tolist():
                neg_q_df = neg_q_df.drop(columns=[mc])
        qc_cols = [col for col in neg_q_df.columns if re_qc.search(col)]
        neg_q_df = neg_q_df[neg_q_df.apply(self.filter_qc, qc_cols=qc_cols, axis=1)]
        rsd_v = self.judge_rsd(neg_q_df, qc_cols)
        if rsd_v < 0.7:
            print('neg定量表qc的rsd小于30的不足0.7，qc用20筛选')
            neg_q_df = neg_q_df[neg_q_df.apply(self.filter_qc, qc_cols=qc_cols, cutoff=20, axis=1)]
            rsd_v = self.judge_rsd(neg_q_df, qc_cols)
            if rsd_v < 0.7:
                if neg_q_df.shape[0] > 20000:
                    print('neg的总峰值超过20000，尝试qc用25筛选')
                    neg_q_df = neg_q_df[neg_q_df.apply(self.filter_qc, qc_cols=qc_cols, cutoff=25, axis=1)]
                    rsd_v = self.judge_rsd(neg_q_df, qc_cols)
                    if rsd_v < 0.7:
                        print_flower_letter('warnning', 'white', 'red', 'anti')
                        print('警告：neg定量表qc用25筛选后的rsd依旧不足0.7  现为%s' % str(rsd_v))
                else:
                    print_flower_letter('warnning', 'white', 'red', 'anti')
                    print('警告：neg定量表qc用20筛选后的rsd依旧不足0.7  现为%s' % str(rsd_v))
            if rsd_v < 0.7:
                for qc in qc_cols:
                    tmp = qc_cols[::]
                    tmp.remove(qc)
                    rsd_tmp = self.judge_rsd(neg_q_df, tmp)
                    if rsd_tmp >= 0.7:
                        print('你可以尝试删除%s' % qc)
                    print('去掉%s后，neg表的rsd为%s' % (qc, rsd_tmp))
        # neg_q_df.insert(0, 'Feature ID', neg_q_df.index + 1)
        print('本项目neg定量表的rsd为%s' % str(rsd_v))

        neg_q_df.insert(0, 'Feature ID', np.arange(0, neg_q_df.shape[0]) + 1)

        # 处理定性文件，增加Theoretical Fragmentation Score ，去重，增加定量信息
        sort_cols = ["Compound", "Compound ID", "Adducts", "Formula", "Fragmentation Score", "Theoretical Fragmentation Score", "Mass Error (ppm)", "Isotope Similarity", "Theoretical Isotope Distribution", "Link", "Description", "Neutral mass (Da)", "m/z", "Charge", "Retention time (min)"]
        quant_col = [0, 1] + list(range(7, pos_q_df.shape[1]))
        with open(self.pos_i, 'rb') as pir:
            try:
                pos_i_df = pd.read_csv(pir, sep=',', quoting=0, quotechar='"', encoding='utf-8').fillna(0)
            except:
                pir.seek(0, 0)
                pos_i_df = pd.read_csv(pir, sep=',', quoting=0, quotechar='"', encoding='GB18030').fillna(0)
        if self.keep_hmdb != 'yes':
            print('不保留HMDB相关的代谢物')
            pos_i_df = pos_i_df[~pos_i_df['Compound ID'].astype(str).str.contain('HMDB')]
        print('正在通过Description筛选定性数据')
        pos_i_df['Description'] = pos_i_df.apply(self.p_object.deal_description,axis=1)
        pos_i_df = pos_i_df[pos_i_df.apply(self.p_object.filter_by_description,axis=1)]
        uniq_ = self.get_unique_index(pos_i_df)
        pos_i_df = pos_i_df.loc[uniq_,]
        # 还需要用description去重
        uniq_ = self.get_unique_index(pos_i_df, use='Description')
        pos_i_df = pos_i_df.loc[uniq_,]
        tmp_s = pos_i_df.apply(self.deal_score, axis=1)
        pos_i_df['Fragmentation Score'] = tmp_s.apply(lambda x: x[0])
        pos_i_df['Theoretical Fragmentation Score'] = tmp_s.apply(lambda x: x[1])

        pos_i_df = pos_i_df[sort_cols]
        pos_i_df = pos_i_df.replace(-np.inf, '_')
        pos_i_df = pd.merge(pos_i_df, pos_q_df.iloc[:, quant_col], on='Compound', how='inner')
        tmp = pos_i_df['Feature ID']
        pos_i_df = pos_i_df.drop(columns='Feature ID')
        pos_i_df.insert(0, 'Feature ID', tmp)

        with open(self.neg_i, 'rb') as nir:
            try:
                neg_i_df = pd.read_csv(nir, sep=',', quoting=0, quotechar='"', encoding='utf-8').fillna(0)
            except:
                nir.seek(0, 0)
                neg_i_df = pd.read_csv(nir, sep=',', quoting=0, quotechar='"', encoding='GB18030').fillna(0)
        if self.keep_hmdb != 'yes':
            neg_i_df = neg_i_df[~neg_i_df['Compound ID'].astype(str).str.contains('HMDB')]
        neg_i_df['Description'] = neg_i_df.apply(self.p_object.deal_description,axis=1)
        neg_i_df = neg_i_df[neg_i_df.apply(self.p_object.filter_by_description,axis=1)]
        uniq_ = self.get_unique_index(neg_i_df)
        neg_i_df = neg_i_df.loc[uniq_,]
        uniq_ = self.get_unique_index(neg_i_df, use='Description')
        neg_i_df = neg_i_df.loc[uniq_,]
        tmp_s = neg_i_df.apply(self.deal_score, axis=1)
        neg_i_df['Fragmentation Score'] = tmp_s.apply(lambda x: x[0])
        neg_i_df['Theoretical Fragmentation Score'] = tmp_s.apply(lambda x: x[1])
        neg_i_df = neg_i_df[sort_cols]
        neg_i_df = neg_i_df.replace(-np.inf, '_')
        quant_col = [0, 1] + list(range(7, neg_q_df.shape[1]))
        neg_i_df = pd.merge(neg_i_df, neg_q_df.iloc[:, quant_col], on='Compound', how='inner')
        tmp = neg_i_df['Feature ID']
        neg_i_df = neg_i_df.drop(columns='Feature ID')
        neg_i_df.insert(0, 'Feature ID', tmp)

        # 导出数据
        # pos_q = self.pos_q + '_dealed.csv'
        # pos_i = self.pos_i + '_dealed.csv'
        # neg_q = self.neg_q + '_dealed.csv'
        # neg_i = self.neg_i + '_dealed.csv'
        pos_q = os.path.join(self.output_dir, os.path.basename(self.pos_q) + "_dealed.csv")
        pos_i = os.path.join(self.output_dir, os.path.basename(self.pos_i) + "_dealed.csv")
        neg_q = os.path.join(self.output_dir, os.path.basename(self.neg_q) + "_dealed.csv")
        neg_i = os.path.join(self.output_dir, os.path.basename(self.neg_i) + "_dealed.csv")
        with open(pos_q, 'w') as pqw:
            pos_q_df.to_csv(pqw, sep=',', quoting=0, quotechar='"', index=False, header=True, encoding='utf-8')
        with open(pos_i, 'wb') as piw:
            tmp_f = StringIO()
            try:
                pos_i_df.to_csv(tmp_f, sep=',', quoting=0, quotechar='"', index=False, header=True, encoding='utf-8')
                piw.write(tmp_f.getvalue().encode('utf8'))
            except:
                pos_i_df.to_csv(tmp_f, sep=',', quoting=0, quotechar='"', index=False, header=True)
                piw.write(tmp_f.getvalue().encode('utf8'))
            tmp_f.close()
        with open(neg_q, 'w') as nqw:
            neg_q_df.to_csv(nqw, sep=',', quoting=0, quotechar='"', index=False, header=True, encoding='utf-8')
        with open(neg_i, 'wb') as niw:
            tmp_f = StringIO()
            try:
                neg_i_df.to_csv(tmp_f, sep=',', quoting=0, quotechar='"', index=False, header=True, encoding='utf-8')
                niw.write(tmp_f.getvalue().encode('utf8'))
            except:
                neg_i_df.to_csv(niw, sep=',', quoting=0, quotechar='"', index=False, header=True)
                niw.write(tmp_f.getvalue().encode('utf8'))
            tmp_f.close()
        # 如果鉴定量小于100要警报，鉴定量是指metlin的大于50pos和neg加起来不到100的
        tmp_pos = pos_i_df['Fragmentation Score'].replace('_', 0)
        tmp_pos = tmp_pos[tmp_pos.astype(float) > 50]
        tmp_neg = neg_i_df['Fragmentation Score'].replace('_', 0)
        tmp_neg = tmp_neg[tmp_neg.astype(float) > 50]
        if tmp_pos.count() + tmp_neg.count() < 100:
            print_flower_letter('warnning', 'white', 'red', 'anti')
            print('鉴定量不足，仅为%s, 请跟生产人员联系' % str(tmp_pos.count() + tmp_neg.count()))

        print('本次pos表峰总量为%s'%str(pos_q_df.shape[0]))
        print('本次neg表峰总量为%s'%str(neg_q_df.shape[0]))

        return pos_i, pos_q, neg_i, neg_q


if __name__ == '__main__':
    # def __init__(self, pos_i, pos_q, neg_i, neg_q, max_row, min_score, max_len, group_table, chname, cutoff_percent)
    parser = argparse.ArgumentParser(description="处理LC代谢的下机数据脚本")
    parser.add_argument("-pos_i", type=str, required=True, help="正离子定性文件")
    parser.add_argument("-pos_q", type=str, required=True, help="正离子定量文件")
    parser.add_argument("-neg_i", type=str, required=True, help="负离子定性文件")
    parser.add_argument("-neg_q", type=str, required=True, help="负离子定量文件")
    parser.add_argument("-filter_free", type=str, default=None, help="不希望筛除的代谢物的Compound列表")
    parser.add_argument("-use_col", type=int, default=0, help="Normalised abundance这一列所在的列数，不填的话就自己查找")


    args = parser.parse_args()
    predeal = PredealLcRaw(args.pos_i,
                           args.pos_q,
                           args.neg_i,
                           args.neg_q,
                           args.use_col,
                           args.filter_free,
                           )
    predeal.run()
