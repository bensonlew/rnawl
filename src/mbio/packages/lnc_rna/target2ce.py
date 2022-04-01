# -*- coding: utf-8 -*-
import os
import numpy as np
import pandas as pd
import argparse
from scipy.stats import pearsonr, spearmanr, kendalltau
from statsmodels.stats.multitest import multipletests
from scipy import stats
import json
from collections import OrderedDict


class Target2ce(object):
    def __init__(self, mrna_target = None,
                 lncrna_target = None,
                 output = None,
                 mirna_family = None,
                 gene_type = None):
        self.mrna_target = mrna_target
        self.lncrna_target = lncrna_target
        self.output = output
        self.mirna_family = mirna_family
        self.ceRNA_list = list()
        self.f2m = dict()
        self.gene_type = gene_type


    def get_ce_mirnas(self):
        target_types = ['mRNA', 'lncRNA']
        target_file = [self.mrna_target, self.lncrna_target]
        ce2mirna = dict()
        ce2type = dict()
        ce2name = dict()
        ce2gene = dict()

        if self.gene_type == "G": 
            cn = 2
        else:
            cn = 1
        for n,target in enumerate(target_file):
            with open(target, 'r') as f:
                f.readline()
                for line in f:
                    cols = line.strip("\n").split("\t")
                    if cols[cn] in ce2mirna:
                        ce2mirna[cols[cn]].add(cols[0])
                    else:
                        ce2mirna[cols[cn]] = set([cols[0]])
                    ce2type[cols[cn]] = target_types[n]
                    ce2name[cols[cn]] = cols[3]
                    if self.gene_type == "G":
                        ce2gene[cols[2]] = cols[2]
                    else:
                        ce2gene[cols[1]] = cols[2]

        return ce2mirna, ce2type, ce2name, ce2gene

    def change_mi2fam(self, ce2m):
        '''
        修改mirna合并family结果, 同时将fam2mirna 字典付给类对象
        '''
        m2f = dict()
        ce2f =dict()
        f2m = dict()
        all_familys = set()
        if self.mirna_family:
            with open(self.mirna_family, 'r') as mf_f:
                mf_f.readline()
                for line in mf_f:
                    cols = line.strip("\n").split("\t")
                    if cols[2] == "-":
                        continue
                    m2f[cols[0]] = cols[2]
                    if cols[2] in f2m:
                        f2m[cols[2]] += "," + cols[0]
                    else:
                        f2m[cols[2]] = cols[0]
                self.f2m = f2m
            for ce, mirnas  in ce2m.items():
                familys = [m2f[x] if x in m2f else x for x in mirnas]
                ce2f[ce] = set(familys)
                all_familys = all_familys | set(familys)
            return ce2f, all_familys
        else:
            for familys in ce2m.values:
                all_familys = all_familys | familys
            return ce2m, all_familys


    def run_ce(self, pvalue_cut):
        ce2mirna, ce2type, ce2name, ce2gene = self.get_ce_mirnas()
        ce2mfamily, all_familys = self.change_mi2fam(ce2mirna)
        ces = ce2mirna.keys()
        all_num = len(all_familys)

        with open(self.output, 'w') as ce_f:
            ce_f.write("ceRNA1\tceRNA2\tceRNA1_name\tceRNA2_name\tceRNA1_type\tceRNA2_type\thit_mirna_num\thit_mirna_familys\tpvalue\tceRNA1_gene\tceRNA2_gene\n")
            for n in range(len(ces) - 1):
                for k in range(n+1, len(ces)):
                    ce1 = ces[n]
                    ce2 = ces[k]
                    ce1_num = len(ce2mfamily[ce1])
                    ce2_num = len(ce2mfamily[ce2])
                    inter_num = len(ce2mfamily[ce1] & ce2mfamily[ce2])
                    if inter_num == 0:
                        continue
                    p_value = stats.hypergeom.pmf(inter_num, all_num, ce1_num, ce2_num)
                    inter_mi_str_list = list()
                    for fam in ce2mfamily[ce1] & ce2mfamily[ce2]:
                        if fam in self.f2m:
                            inter_mi_str = fam + ":" + self.f2m[fam]
                        else:
                            inter_mi_str = fam
                        inter_mi_str_list.append(inter_mi_str)

                    if p_value < pvalue_cut and str(inter_num)>=3:
                        ce_f.write("\t".join([
                            ce1, ce2, ce2name[ce1], ce2name[ce2], ce2type[ce1], ce2type[ce2],
                            str(inter_num), ";".join(inter_mi_str_list),  str(p_value),
                            ce2gene[ce1], ce2gene[ce2]
                        ]) + "\n")

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description='for gene correlation calculation and adjustment')
    parser.add_argument('-mrna_target', type=str, default=None,
                        help="mirna target matrix.")
    parser.add_argument('-lncrna_target', type=str, default=None,
                        help="mirna target matrix.")

    parser.add_argument('-output', type=str, default="ce.xls", help='output directory.')
    parser.add_argument('-gene_type', type=str, default="G", help='ce target is gene or transcript')
    parser.add_argument('-mirna_family', type=str, default=None, help='miRNA family file')
    parser.add_argument('-pvalue_cutoff', type=float, default=0.05, help='pvalue cutoff. Default: 0.05')
    parser.add_argument('-anno', type=str, default=None, help='the dir of annotation')
    # ----------------------------------------------------------------------------------------------
    args = parser.parse_args()
    toolbox = Target2ce( mrna_target = args.mrna_target,
                         lncrna_target = args.lncrna_target,
                         output = args.output,
                         mirna_family = args.mirna_family,
                         gene_type = args.gene_type
    )
    toolbox.run_ce(pvalue_cut=args.pvalue_cutoff)

