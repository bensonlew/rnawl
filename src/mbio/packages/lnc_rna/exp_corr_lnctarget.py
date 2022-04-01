# -*- coding: utf-8 -*-
import os
import numpy as np
import pandas as pd
import argparse
from scipy.stats import pearsonr, spearmanr, kendalltau
from statsmodels.stats.multitest import multipletests
import json
from collections import OrderedDict
from multiprocessing import Pool

df2 = pd.DataFrame()
global cor_cutoff
cor_cutoff = 0.8



def sub_pearson(lnc_exp):
    global df2
    global cor_cutoff
    lncrna_id = lnc_exp[0]
    lnc_exp1 = lnc_exp[1]
    sub_corr = list()
    for j in df2.index:
        target_id = j
        try:
            target_exp = df2.ix[target_id]
            corrtest = list(pearsonr(list(lnc_exp1), list(target_exp)))
            if np.isnan(corrtest[1]):
                corrtest[1] = 1.0
            # print corrtest
            # print corrtest
            if not np.isnan(corrtest[0]):
                if abs(corrtest[0]) > abs(cor_cutoff):
                    sub_corr.append([lncrna_id, target_id, corrtest[0], corrtest[1]])
        except Exception,e:
            print e
    return sub_corr

def sub_spearman(lnc_exp):
    global df2
    global cor_cutoff
    lncrna_id = lnc_exp[0]
    lnc_exp1 = lnc_exp[1]
    sub_corr = list()
    for j in df2.index:
        target_id = j
        try:
            target_exp = df2.ix[target_id]
            corrtest = list(spearmanr(list(lnc_exp[1]), list(target_exp)))
            # print corrtest
            if np.isnan(corrtest[1]):
                corrtest[1] = 1.0
            
            if not np.isnan(corrtest[0]):
                if abs(corrtest[0]) > abs(cor_cutoff):
                    sub_corr.append([lncrna_id, target_id, corrtest[0], corrtest[1]])
        except Exception,e:
            print e
    return sub_corr

def sub_kendall(lnc_exp):
    global df2
    global cor_cutoff
    lncrna_id = lnc_exp[0]
    lnc_exp1 = lnc_exp[1]
    sub_corr = list()
    for j in df2.index:
        target_id = j
        try:
            target_exp = df2.ix[target_id]
            corrtest = list(kendalltau(list(lnc_exp[1]), list(target_exp)))
            # print corrtest
            if np.isnan(corrtest[1]):
                corrtest[1] = 1.0
            if not np.isnan(corrtest[0]):
                if abs(corrtest[0]) > abs(cor_cutoff):
                    sub_corr.append([lncrna_id, target_id, corrtest[0], corrtest[1]])
        except Exception,e:
            print e
    return sub_corr



class Expresscor(object):
    def __init__(self, padjust_way="fdr_bh", corr_way="pearson",
                 output="corr.xls", exp_matrix=None, exp_target=None, corr_cutoff=0,
                 pvalue_cutoff=None, qvalue_cutoff=None, sig_type=None,
                 g_or_t=None ,anno=None):
        self.exp = exp_matrix
        self.target_exp = exp_target
        self.corr_cutoff = corr_cutoff
        self.qvalue_cutoff = qvalue_cutoff
        self.pvalue_cutoff = pvalue_cutoff
        if padjust_way == "bh":
            self.padjust_way = "fdr_bh"
        elif padjust_way == "by":
            self.padjust_way = "fdr_by"
        else:
            self.padjust_way = padjust_way
        self.corr_way = corr_way
        global cor_cutoff
        cor_cutoff = self.corr_cutoff
        '''
        if self.corr_way == "pearson":
            self.corr_way = pearsonr
        if self.corr_way == "spearman":
            self.corr_way = spearmanr
        if self.corr_way == "spearmanr":
            self.corr_way = spearmanr
        if self.corr_way == "kendall":
            self.corr_way = kendalltau
        '''
        self.output = output
        self.sig_type = sig_type
        self.g_or_t = g_or_t
        self.anno = anno

    def calculate_corr(self, pvalue_cutoff=None, qvalue_cutoff=None, thread=20):
        exp_matrix = self.exp
        df1 = pd.read_table(exp_matrix, sep="\t", header=0, index_col=0)
        df1 = df1[(df1 != 0).any(axis=1)]
        if "is_new" in df1.columns:
            df1 = df1.drop(columns=["is_new"])
        if "rna_type" in df1.columns:
            df1 = df1.drop(columns=["rna_type"])
        del df1.index.name
        # df1 = df1.transpose()

        df2 = pd.read_table(self.target_exp, sep="\t", header=0, index_col=0)
        df2 = df2[(df2 != 0).any(axis=1)]
        if "is_new" in df2.columns:
            df2 = df2.drop(columns=["is_new"])
        if "rna_type" in df2.columns:
            df2 = df2.drop(columns=["rna_type"])
        del df2.index.name

        # 过滤矩阵， 过滤低表达的和偏差小的数据

        if len(df1.columns) == 2:
            max_cut = 10
            stc_cut = 0.25
        else:
            max_cut = 1
            stc_cut = 0.1
        df1_max = df1.max(axis=1)
        df1_max_set = set(df1_max[df1_max>max_cut].index)
        df1_std = df1.std(axis=1)
        df1_mean = df1.mean(axis=1)
        df1_stc_set = set(df1[df1_std/df1_mean > stc_cut].index)
        df1 = df1.loc[df1_max_set.intersection(df1_stc_set), :]
        df2_max = df2.max(axis=1)
        df2_max_set = set(df2_max[df2_max>max_cut].index)
        df2_std = df2.std(axis=1)
        df2_mean = df2.mean(axis=1)
        df2_stc_set = set(df2[df2_std/df2_mean > stc_cut].index)
        df2 = df2.loc[df2_max_set.intersection(df2_stc_set), :]


        # df1.to_csv(self.output + "test.lnc", sep="\t", index=False)
        # df2.to_csv(self.output + "test.m", sep="\t", index=False)
        global df2

        lnc_exp_list = list()
        for i in df1.index:
            lncrna_id = i
            lnc_exp = list(df1.ix[lncrna_id])
            lnc_exp_list.append([lncrna_id, lnc_exp])

        p = Pool(thread)
        # print lnc_exp_list

        if self.corr_way == "pearson":
            corr_list = p.map(sub_pearson, lnc_exp_list)
        elif self.corr_way == "kendall":
            corr_list = p.map(sub_kendall, lnc_exp_list)
        else:
            corr_list = p.map(sub_spearman, lnc_exp_list)

        p.close()
        p.join()
        corr = list()
        for sub_corr in corr_list:
            corr.extend(sub_corr)

        if len(corr) > 0:
            corr_pd = pd.DataFrame(corr)
            corr_pd.rename(columns = { 0:"lncrna", 1: "target_gene", 2: "corr", 3: "pvalue"}, inplace=True)
            # corr_pd.to_csv("corr.xls", sep="\t", index=False)
        else:
            corr_pd = pd.DataFrame({"lnc_rna":[], "target_gene": [], "corr": [], "pvalue": []})

        # corr_pd.to_csv(self.output + "test.corr", sep="\t", index=False)
        if len(corr) > 0:
            flat_adjust_pvalue = multipletests(corr_pd["pvalue"], alpha=1, method=self.padjust_way, is_sorted=False, returnsorted=False)
            corr_pd['qvalue'] = flat_adjust_pvalue[1]
        else:
            corr_pd['qvalue'] = [1 for p in corr_pd["pvalue"]]

        if self.pvalue_cutoff:
            corr_filter = corr_pd[(corr_pd['pvalue']<=self.pvalue_cutoff)]
        else:
            pass

        if self.qvalue_cutoff:
            corr_filter = corr_pd[(corr_pd['qvalue']<=self.qvalue_cutoff)]
        else:
            pass

        corr_filter.to_csv(self.output, sep="\t", index=False)


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description='for gene correlation calculation and adjustment')
    parser.add_argument('-exp', type=str, default=None,
                        help="mirna express matrix.")
    parser.add_argument('-exp_target', type=str, default=None,
                        help="gene or transcript express matrix.")
    parser.add_argument('-output', type=str, default="corr.xls", help='output directory.')
    parser.add_argument('-pvalue_cutoff', type=float, default=0.05, help='pvalue cutoff. Default: 0.05')
    parser.add_argument('-qvalue_cutoff', type=float, default=0.05, help='qvalue cutoff. Default: 0.05')
    parser.add_argument('-cor_cutoff', type=float, default=0, help='correlation cutoff. Default: 0.05')
    parser.add_argument('-corr_way', type=str, default="pearson")
    parser.add_argument('-g_or_t', type=str, default=None, help='the gene level or transcript level')
    parser.add_argument('-anno', type=str, default=None, help='the dir of annotation')
    parser.add_argument('-sig_type', type=int, default=1, help='1 means we '
                                                               'use qvalue '
                                                               'to select, '
                                                               'while 0 '
                                                               'means pvalue to select')
    parser.add_argument('-padjust_way', type=str, default="fdr_bh",
                        help='http://www.statsmodels.org/devel/generated/statsmodels.stats.multitest.multipletests.html'
                             'bonferroni, fdr_bh : Benjamini/Hochberg (non-negative), fdr_by : Benjamini/Yekutieli (negative)')

    # ----------------------------------------------------------------------------------------------
    args = parser.parse_args()
    toolbox = Expresscor(exp_matrix=args.exp,
                         exp_target=args.exp_target,
                         output=args.output,
                         sig_type=args.sig_type,
                         corr_cutoff=args.cor_cutoff,
                         pvalue_cutoff=args.pvalue_cutoff,
                         qvalue_cutoff=args.qvalue_cutoff,
                         padjust_way=args.padjust_way,
                         corr_way=args.corr_way,
                         g_or_t=args.g_or_t,
                         anno=args.anno)
    toolbox.calculate_corr(pvalue_cutoff=args.pvalue_cutoff, qvalue_cutoff=args.qvalue_cutoff)

