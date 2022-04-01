# -*- coding: utf-8 -*-
import os
import numpy as np
import pandas as pd
import argparse
from scipy.stats import pearsonr, spearmanr, kendalltau
from statsmodels.stats.multitest import multipletests
import json
from collections import OrderedDict


class ExpressCe(object):
    def __init__(self, padjust_way="fdr_bh", corr_way="pearson",
                 output="corr.xls",
                 exp_mirna=None,
                 exp_cerna=None,
                 corr_mirna2cerna=None,
                 ce_pair=None,
                 cor_cutoff=0,
                 pvalue_cutoff=None, qvalue_cutoff=None, sig_type=None,
                 g_or_t=None ,anno=None):
        self.exp_mirna = exp_mirna
        self.exp_cerna = exp_cerna
        self.ce_pair = ce_pair
        self.corr_mirna2cerna = corr_mirna2cerna
        self.cor_cutoff = cor_cutoff
        self.qvalue_cutoff = qvalue_cutoff
        self.pvalue_cutoff = pvalue_cutoff
        if padjust_way == "bh":
            self.padjust_way = "fdr_bh"
        elif padjust_way == "by":
            self.padjust_way = "fdr_by"
        else:
            self.padjust_way = padjust_way
        self.corr_way = corr_way
        if self.corr_way == "pearson":
            self.corr_way = pearsonr
        if self.corr_way == "spearman":
            self.corr_way = spearmanr
        if self.corr_way == "kendall":
            self.corr_way = kendalltau
        self.output = output
        self.sig_type = sig_type
        self.g_or_t = g_or_t
        self.anno = anno

    def get_cerna2mirna_corr(self):
        self.mice2corr = dict()
        if self.corr_mirna2cerna:
            for f_path in self.corr_mirna2cerna.split(","):
                with open(f_path, 'r') as f:
                    file_corr = False
                    head = f.readline()
                    if head.split("\t")[-3] == "corr":
                        file_corr = True
                    for line in f:
                        cols = line.strip().split("\t")
                        if self.g_or_t == "G":
                            cn = 2
                        else:
                            cn = 1
                        if file_corr:
                            self.mice2corr[cols[0] + "__" + cols[cn]] = cols[-3]
                        else:
                            self.mice2corr[cols[0] + "__" + cols[cn]] = "0"


    def add_mirna_corr(self, mirna, ce1_id, ce2_id):
        '''
        mirna添加相关系数
        过滤没有相关系数的miRNA对
        '''
        pair_str = list()
        regu1 = False
        regu2 = False
        corr_list = list()
        for fam in mirna.split(";"):
            mi_ce_corrs = list()
            for mirna in fam.split(":")[-1].split(","):
                if mirna + "__" + ce1_id in self.mice2corr:
                    cor = self.mice2corr[mirna + "__" + ce1_id]
                    if float(cor) <= 0:
                        regu1 = True
                else:
                    continue
                    # cor = "0"
                if mirna + "__" + ce2_id in self.mice2corr:
                    cor1 = cor
                    cor2 = self.mice2corr[mirna + "__" + ce2_id]
                    cor = cor + "|" + cor2
                    if float(self.mice2corr[mirna + "__" + ce2_id]) <= 0:
                        regu2 = True
                        corr_list.append((float(cor1), float(cor2)))
                else:
                    continue
                    # cor = cor + "|" + "0"
                mi_ce_corrs.append(mirna + "(" + cor + ")")
            if len(mi_ce_corrs) != 0:
                if len(fam.split(":")) == 2:
                    mi_ce_str = fam.split(":")[-2] + ":" + ",".join(mi_ce_corrs)
                else:
                    mi_ce_str = ",".join(mi_ce_corrs)
                pair_str.append(mi_ce_str)

        return pair_str, regu1 & regu2, corr_list

    def calculate_corr(self, pvalue_cutoff=None, qvalue_cutoff=None):
        self.get_cerna2mirna_corr()
        exp_matrix = self.exp_cerna
        ce_exp = pd.read_table(exp_matrix, sep="\t", header=0, index_col=0)
        '''
        df1 = df1[(df1 != 0).any(axis=1)]
        del df1.index.name
        # df1 = df1.transpose()
        df2 = pd.read_table(self.target_exp, sep="\t", header=0, index_col=0)
        df2 = df2[(df2 != 0).any(axis=1)]
        
        del df2.index.name
        '''
        # df2 = df2.transpose()
        df3 = pd.read_table(self.ce_pair, sep="\t", header=0)
        coeffmat = np.zeros(df3.shape[0])
        pvalmat = np.zeros(df3.shape[0])
        ce_mirnas = list()
        is_regus = list()
        regu_similar = list()
        sens_corr = list()
        mi_num = list()
        for i in range(df3.shape[0]):
            ce1_id = df3.loc[i][0]
            ce2_id = df3.loc[i][1]

            try:
                # print df1.ix[mirna_id]
                # print df2.ix[gene_id]
                if len(ce_exp.ix[ce1_id].shape) == 2:
                    ce1_exp = ce_exp.ix[ce1_id].iloc[0]
                else:
                    ce1_exp = ce_exp.ix[ce1_id]

                if len(ce_exp.ix[ce2_id].shape) == 2:
                    ce2_exp = ce_exp.ix[ce2_id].iloc[0]
                else:
                    ce2_exp = ce_exp.ix[ce2_id]

                corrtest = self.corr_way(list(ce1_exp), list(ce2_exp))

                if np.isnan(corrtest[0]):
                    coeffmat[i] = 0
                else:
                    coeffmat[i] = corrtest[0]
                if np.isnan(corrtest[1]):
                    pvalmat[i] = 1
                else:
                    pvalmat[i] = corrtest[1]
            except:
                coeffmat[i] = 0
                pvalmat[i] = 1

            mirna = df3.loc[i]['hit_mirna_familys']

            pair_str, is_regu, corr_list = self.add_mirna_corr(mirna, ce1_id, ce2_id)
            #计算 Regulation similarity 和 Sensitivity correlation
            print "corr_list is {}".format(corr_list)
            if(len(corr_list) != 0):
                m = len(corr_list)
                ce_cor = coeffmat[i]
                sim = 1 - sum([abs((a-b)/(a+b))**m for a,b in corr_list if a+b != 0])/m
                try:
                    corr = ce_cor - sum([(ce_cor - (a*b))/(((1-a*a)**0.5)*((1-b*b)**0.5)) for a,b in corr_list])/m
                except:
                    corr = ce_cor
            else:
                sim = 0
                corr = 0
                m = 0

            regu_similar.append(sim)
            sens_corr.append(corr)

            ce_mirnas.append(";".join(pair_str))
            mi_num.append(len(pair_str))
            is_regus.append(is_regu)

        if len(pvalmat) != 0:
            flat_adjust_pvalue = multipletests(pvalmat, alpha=1, method=self.padjust_way, is_sorted=False, returnsorted=False)
            df3['corr'] = coeffmat
            df3['corr_pvalue'] = pvalmat
            df3['corr_qvalue'] = flat_adjust_pvalue[1]
            df3['hit_mirna_corr'] = ce_mirnas
            df3['is_regu'] = is_regus
            df3['regu_similar'] = regu_similar
            df3['sens_corr'] = sens_corr
        else:
            df3['corr'] = []
            df3['corr_pvalue'] = pvalmat
            df3['corr_qvalue'] = []
            df3['hit_mirna_corr'] = ce_mirnas
            df3['is_regu'] = is_regus
            df3['regu_similar'] = regu_similar
            df3['sens_corr'] = sens_corr
        df3['hit_mirna_num'] = mi_num

        df3.to_csv(self.output + ".unfilter", sep="\t", index=False)

        df3_filter = df3[df3['corr'] > self.cor_cutoff]
        if self.corr_mirna2cerna:
            df3_filter = df3_filter[df3_filter['is_regu'] == True]

        if self.pvalue_cutoff:
            df3_filter = df3_filter[(df3_filter['corr_pvalue']<self.pvalue_cutoff)]
        else:
            pass

        if self.qvalue_cutoff:
            df3_filter = df3_filter[(df3_filter['corr_qvalue']<self.qvalue_cutoff)]
        else:
            pass

        df_filter = df3_filter[['ceRNA1', 'ceRNA2', 'ceRNA1_name', 'ceRNA2_name', 'ceRNA1_type', 'ceRNA2_type', 'hit_mirna_num', 'pvalue', 'corr', 'corr_pvalue', 'corr_qvalue', 'hit_mirna_corr', 'ceRNA1_gene', 'ceRNA2_gene', 'regu_similar', 'sens_corr']]
        df_filter.to_csv(self.output, sep="\t", index=False)


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description='for gene correlation calculation and adjustment')
    parser.add_argument('-exp_mirna', type=str, default=None,
                        help="mirna express matrix.")
    parser.add_argument('-exp_cerna', type=str, default=None,
                        help="gene or transcript express matrix.")
    parser.add_argument('-corr_mirna2cerna', type=str, default=None,
                        help="corr matrix of mirna2cerna.")
    parser.add_argument('-ce_pair', type=str, default=None,
                        help="ce rna predict table.")
    parser.add_argument('-output', type=str, default="corr.xls", help='output directory.')
    parser.add_argument('-pvalue_cutoff', type=float, default=1, help='pvalue cutoff. Default: 1')
    parser.add_argument('-qvalue_cutoff', type=float, default=1, help='qvalue cutoff. Default: 1')
    parser.add_argument('-cor_cutoff', type=float, default=0, help='correlation cutoff. Default: 0')
    parser.add_argument('-corr_way', type=str, default="pearson")
    parser.add_argument('-g_or_t', type=str, default="G", help='the gene level or transcript level')
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
    toolbox = ExpressCe(exp_mirna=args.exp_mirna,
                        exp_cerna=args.exp_cerna,
                        corr_mirna2cerna=args.corr_mirna2cerna,
                        ce_pair=args.ce_pair,
                        output=args.output,
                        sig_type=args.sig_type,
                        cor_cutoff=args.cor_cutoff,
                        pvalue_cutoff=args.pvalue_cutoff,
                        qvalue_cutoff=args.qvalue_cutoff,
                        padjust_way=args.padjust_way,
                        corr_way=args.corr_way,
                        g_or_t=args.g_or_t,
                        anno=args.anno)
    toolbox.calculate_corr(pvalue_cutoff=args.pvalue_cutoff, qvalue_cutoff=args.qvalue_cutoff)

