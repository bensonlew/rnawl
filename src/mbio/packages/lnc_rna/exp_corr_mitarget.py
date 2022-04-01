# -*- coding: utf-8 -*-
import os
import numpy as np
import pandas as pd
import argparse
from scipy.stats import pearsonr, spearmanr, kendalltau
from statsmodels.stats.multitest import multipletests
import json
from collections import OrderedDict


class Expresscor(object):
    def __init__(self, padjust_way="fdr_bh", corr_way="pearson",
                 output="corr.xls", exp_matrix=None, exp_target=None, rna_target=None, cor_cutoff=0,
                 pvalue_cutoff=None, qvalue_cutoff=None, sig_type=None,
                 g_or_t=None ,anno=None):
        self.exp = exp_matrix
        self.target_exp = exp_target
        self.rna_target = rna_target
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

    def calculate_corr(self, pvalue_cutoff=None, qvalue_cutoff=None):
        exp_matrix = self.exp
        df1 = pd.read_table(exp_matrix, sep="\t", header=0, index_col=0)
        df1 = df1[(df1 != 0).any(axis=1)]
        del df1.index.name
        # df1 = df1.transpose()
        df2 = pd.read_table(self.target_exp, sep="\t", header=0, index_col=0)
        df2 = df2[(df2 != 0).any(axis=1)]
        del df2.index.name
        # df2 = df2.transpose()
        df3 = pd.read_table(self.rna_target, sep="\t", header=0)
        coeffmat = np.zeros(df3.shape[0])
        pvalmat = np.zeros(df3.shape[0])
        for i in range(df3.shape[0]):
            mirna_id = df3.loc[i][0]
            if self.g_or_t == "G":
                gene_id = df3.loc[i][2]
            else:
                gene_id = df3.loc[i][1]
            
            try:
                # print df1.ix[mirna_id]
                # print df2.ix[gene_id]
                if len(df1.ix[mirna_id].shape) == 2:
                    mi_exp = df1.ix[mirna_id].iloc[0]
                else:
                    mi_exp = df1.ix[mirna_id]

                # print "***"
                # print mi_exp
                # print df2.ix[gene_id]
                #print corrtest
                #print list(mi_exp), list(df2.ix[gene_id])
                corrtest = self.corr_way(list(mi_exp), list(df2.ix[gene_id]))
                #print corrtest
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
        # print pvalmat
        if len(pvalmat) != 0:
            flat_adjust_pvalue = multipletests(pvalmat, alpha=1, method=self.padjust_way, is_sorted=False, returnsorted=False)
            df3['corr'] = coeffmat
            df3['pvalue'] = pvalmat
            df3['qvalue'] = flat_adjust_pvalue[1]
        else:
            df3['corr'] = []
            df3['pvalue'] = []
            df3['qvalue'] = []


        df3_filter = df3[df3['corr'] < - abs(self.cor_cutoff)]

        # df3_filter = df3

        if self.pvalue_cutoff:
            df3_filter = df3_filter[(df3_filter['pvalue']<=self.pvalue_cutoff)]
        else:
            pass

        if self.qvalue_cutoff:
            df3_filter = df3_filter[(df3_filter['qvalue']<=self.qvalue_cutoff)]
        else:
            pass

        df_filter = df3_filter[['query', 'target', 'gene_id', 'gene_name', 'gene_description', 'mirtarbase', 'score', 'energy', 'corr', 'pvalue', 'qvalue']]

        df_filter.to_csv(self.output, sep="\t", index=False)


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description='for gene correlation calculation and adjustment')
    parser.add_argument('-exp', type=str, default=None,
                        help="mirna express matrix.")
    parser.add_argument('-exp_target', type=str, default=None,
                        help="gene or transcript express matrix.")
    parser.add_argument('-rna_target', type=str, default=None,
                        help="mirna target matrix.")
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
                         rna_target=args.rna_target,
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

