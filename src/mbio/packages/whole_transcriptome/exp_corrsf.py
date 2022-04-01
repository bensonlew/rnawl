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
    def __init__(self, padjust_way="fdr_bh", corr_way="spearmanr",
                 output=None, exp_matrix=None, cor_cutoff=None,category=None,
                 pvalue_cutoff=None, qvalue_cutoff=None, sig_type=None,
                 g_or_t=None ,anno=None):
        self.exp = exp_matrix
        self.cor_cutoff = cor_cutoff
        self.qvalue_cutoff = qvalue_cutoff
        self.pvalue_cutoff = pvalue_cutoff
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
        self.category=category

    def calculate_corr(self, pvalue_cutoff=None, qvalue_cutoff=None):
        exp_matrix = self.exp
        df1 = pd.read_table(exp_matrix, sep="\t", header=0, index_col=0)
        df1 = df1[(df1 != 0).any(axis=1)]
        del df1.index.name
        df1 = df1.transpose()
        df2 = pd.read_table(exp_matrix, sep="\t", header=0, index_col=0)
        df2 = df2[(df2 != 0).any(axis=1)]
        del df2.index.name
        df2 = df2.transpose()
        coeffmat = np.zeros((df1.shape[1], df2.shape[1]))
        pvalmat = np.zeros((df1.shape[1], df2.shape[1]))
        for i in range(df1.shape[1]):
            for j in range(df2.shape[1]):
                corrtest = self.corr_way(df1[df1.columns[i]], df2[df2.columns[j]])
                coeffmat[i,j] = corrtest[0]
                pvalmat[i,j] = corrtest[1]
        dfcoeff = pd.DataFrame(coeffmat, columns=df2.columns, index=df1.columns)
        dfpvals = pd.DataFrame(pvalmat, columns=df2.columns, index=df1.columns)
        flat_pvalue = [item for sublist in dfpvals.values.tolist() for item in sublist]
        flat_adjust_pvalue = multipletests(flat_pvalue, alpha=0.05, method=self.padjust_way, is_sorted=False, returnsorted=False)
        adjust_p_list = flat_adjust_pvalue[1].tolist()
        size = df1.shape[1]
        nested_adjust_pvalue = [adjust_p_list[i:i+size] for i  in range(0, len(adjust_p_list), size)]
        df_adjust =  pd.DataFrame(nested_adjust_pvalue, index=df1.columns, columns=df1.columns)
        output = self.output
        if output is not None:
            if not os.path.exists(output):
                os.mkdir(output)
        else:
            output = os.getcwd()
        out_dir_file = os.path.join(output, "express_correlation_info_pre")
        with open(out_dir_file, "w") as f4:
            seq = ("gene_id_1", "gene_id_2", "cor", "p_value", "q_value")
            f4.write("\t".join(seq) + "\n")
            for i in range(df1.shape[1]):
                for j in range(df1.shape[1]):
                    seq1 = [dfcoeff.index[i], dfcoeff.columns[j], str(dfcoeff.loc[dfcoeff.index[i], dfcoeff.columns[j]])]
                    seq2 = str(dfpvals.loc[dfpvals.index[i], dfpvals.columns[j]])
                    seq3 = str(df_adjust.loc[df_adjust.index[i], df_adjust.columns[j]])
                    f4.write("\t".join(seq1) + "\t" + seq2 + "\t" + seq3 + "\n")
        df_info = pd.read_table(out_dir_file, header=0, sep="\t")
        out_dir_file2 = os.path.join(output, "express_correlation_info_temp.xls")
        out_dir_file3 = os.path.join(output, "express_correlation_info.xls")
        if not self.sig_type:
            df_selected_temp = df_info[(abs(df_info["cor"]) > self.cor_cutoff) & (df_info["p_value"] < self.pvalue_cutoff)]
            df_selected_temp1 = df_selected_temp[df_selected_temp['gene_id_1'] != df_selected_temp['gene_id_2']]
            df_selected_temp1.to_csv(out_dir_file2, sep="\t", index=False)
            with open(out_dir_file2) as f_tmp, open(out_dir_file3, "w") as f_tmp1:
                temp_list = list()
                f_tmp1.write(f_tmp.readline())
                for line in f_tmp:
                    a, b = line.strip().split("\t")[0:2]
                    if [a, b] not in temp_list:
                        temp_list.append([a, b])
                        if not [b, a] in temp_list:
                            f_tmp1.write(line)
            df_selected = pd.read_table(out_dir_file3, header=0, sep="\t")
        if self.sig_type:
            df_selected_temp = df_info[(abs(df_info["cor"]) > self.cor_cutoff) & (df_info["q_value"] < self.qvalue_cutoff)]
            df_selected_temp1 = df_selected_temp[df_selected_temp['gene_id_1'] != df_selected_temp['gene_id_2']]
            df_selected_temp1.to_csv(out_dir_file2, sep="\t", index=False)
            with open(out_dir_file2) as f_tmp, open(out_dir_file3, "w") as f_tmp1:
                temp_list = list()
                f_tmp1.write(f_tmp.readline())
                for line in f_tmp:
                    a, b = line.strip().split("\t")[0:2]
                    if [a, b] not in temp_list:
                        temp_list.append([a, b])
                        if not [b, a] in temp_list:
                            f_tmp1.write(line)
            df_selected = pd.read_table(out_dir_file3, header=0, sep="\t")
        # package的print会显示在.o文件里面
        # anno_xls = os.path.join(self.anno, 'refannot_class/all_annot.xls')
        df_anno = pd.read_table(self.anno)
        if self.g_or_t.lower() == "g":
            df_gene_pre = df_anno.loc[:, ['gene_id', 'gene_name']]
            idx = df_gene_pre[df_gene_pre['gene_name'].isnull()].index
            df_gene_pre.loc[idx, 'gene_name'] = df_gene_pre.loc[idx, 'gene_id']
            df_gene_tmp = df_gene_pre.drop_duplicates('gene_id', keep='first', inplace=False)
            id2name_dict = OrderedDict(zip(df_gene_tmp .gene_id, df_gene_tmp .gene_name))
            df_source_pre = df_selected.loc[:, ['gene_id_1', 'gene_id_2', 'cor']]
            df_source_pre.rename(columns={'gene_id_1': 'id', 'gene_id_2': 'name'}, inplace=True)
            df_nodes = df_source_pre[['id', 'name']]
            list1 = df_nodes['id'].tolist()
            list2 = df_nodes['name'].tolist()
            whole_list_dict = OrderedDict()
            whole_list = list(OrderedDict.fromkeys(list1 + list2))
            whole_name_list = [id2name_dict.get(x, x) for x in whole_list]
            for i in range(len(whole_list)):
                whole_list_dict[whole_list[i]] = i
            df_source_column1 =  [whole_list_dict[x] for x in list1]
            df_source_column2 = [whole_list_dict[x] for x in list2]
            df_source_column3 = df_source_pre['cor']
            df_links = pd.DataFrame(OrderedDict({'source': df_source_column1, 'target': df_source_column2, 'distance': df_source_column3}))
            new_list = list()
            # nodes_links_dict['links'] = df_links.to_dict('records')
            for i in df_links.to_dict('records'):
                new_dict = OrderedDict()
                for key in i.keys():
                    if key != "distance":
                        value = int(i[key])
                        new_dict.update({key: value})
                    else:
                        new_dict.update({"distance": i['distance']})
                new_list.append(new_dict)
                # new_dict.clear()
                del new_dict
            id2type = dict()
            with open(self.category) as ca:
                for line in ca.readlines():
                    line = line.strip().split("\t")
                    id2type[line[0]] = line[1]
            nodes_links_dict = OrderedDict()
            group = [id2type[i] for i in whole_list]
            df_nodes_new_pre =  pd.DataFrame(OrderedDict({"id": whole_list, "name": whole_name_list, "group": group}))
            nodes_links_dict['nodes'] = df_nodes_new_pre.to_dict('records')
            nodes_links_dict['links'] = new_list
            out_dir_file3 = os.path.join(output, "record.json")
            with open(out_dir_file3, "w") as f:
                json.dump(nodes_links_dict, f, separators=(',', ': '), indent=2)
        if self.g_or_t.lower() == "t":
            df_gene_pre = df_anno.loc[:, ['transcript_id', 'gene_name']]
            idx = df_gene_pre[df_gene_pre['gene_name'].isnull()].index
            df_gene_pre.loc[idx, 'gene_name'] = df_gene_pre.loc[idx, 'transcript_id']
            df_gene_tmp = df_gene_pre.drop_duplicates('transcript_id', keep='first', inplace=False)
            id2name_dict = OrderedDict(zip(df_gene_tmp.transcript_id, df_gene_tmp.gene_name))
            df_selected_test = df_selected
            df_selected_test.rename(columns={'gene_id_1': 'transcript_id_1', 'gene_id_2': 'transcript_id_2'}, inplace=True)
            df_source_pre = df_selected_test.loc[:, ['transcript_id_1', 'transcript_id_2', 'cor']]
            df_source_pre.rename(columns={'transcript_id_1': 'id', 'transcript_id_2': 'name'}, inplace=True)
            df_nodes = df_source_pre[['id', 'name']]
            list1 = df_nodes['id'].tolist()
            list2 = df_nodes['name'].tolist()
            whole_list_dict = OrderedDict()
            whole_list = list(OrderedDict.fromkeys(list1 + list2))
            whole_name_list = [id2name_dict[x] for x in whole_list]
            for i in range(len(whole_list)):
                whole_list_dict[whole_list[i]] = i
            df_source_column1 = [whole_list_dict[x] for x in list1]
            df_source_column2 = [whole_list_dict[x] for x in list2]
            df_source_column3 = df_source_pre['cor']
            df_links = pd.DataFrame(
                OrderedDict({'source': df_source_column1, 'target': df_source_column2, 'distance': df_source_column3}))
            new_list = list()
            # nodes_links_dict['links'] = df_links.to_dict('records')
            for i in df_links.to_dict('records'):
                new_dict = OrderedDict()
                for key in i.keys():
                    if key != "distance":
                        value = int(i[key])
                        new_dict.update({key: value})
                    else:
                        new_dict.update({"distance": i['distance']})
                new_list.append(new_dict)
                # new_dict.clear()
                del new_dict
            nodes_links_dict = OrderedDict()
            id2type = dict()
            with open(self.category) as ca:
                for line in ca.readlines():
                    line = line.strip().split("\t")
                    id2type[line[0]] = line[1]
            nodes_links_dict = OrderedDict()
            group = [id2type[i] for i in whole_list]
            df_nodes_new_pre = pd.DataFrame(OrderedDict({"id": whole_list, "name": whole_name_list, "group": group}))
            nodes_links_dict['nodes'] = df_nodes_new_pre.to_dict('records')
            nodes_links_dict['links'] = new_list
            out_dir_file3 = os.path.join(output, "record.json")
            with open(out_dir_file3, "w") as f:
                json.dump(nodes_links_dict, f, separators=(',', ': '), indent=2)

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description='for gene correlation calculation and adjustment')
    parser.add_argument('-exp', type=str, default=None,
                        help="gene or transcript express matrix.")
    parser.add_argument('-output', type=str, default=None, help='output directory.')
    parser.add_argument('-category', type=str, default=None, help='id to category')
    parser.add_argument('-pvalue_cutoff', type=float, default=0.05, help='pvalue cutoff. Default: 0.05')
    parser.add_argument('-qvalue_cutoff', type=float, default=0.05, help='qvalue cutoff. Default: 0.05')
    parser.add_argument('-cor_cutoff', type=float, default=0.8, help='correlation cutoff. Default: 0.05')
    parser.add_argument('-corr_way', type=str, default="spearmanr")
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
                         output=args.output,
                         sig_type=args.sig_type,
                         cor_cutoff=args.cor_cutoff,
                         pvalue_cutoff=args.pvalue_cutoff,
                         qvalue_cutoff=args.qvalue_cutoff,
                         padjust_way=args.padjust_way,
                         corr_way=args.corr_way,
                         g_or_t=args.g_or_t,
                         category=args.category,
                         anno=args.anno)
    toolbox.calculate_corr(pvalue_cutoff=args.pvalue_cutoff, qvalue_cutoff=args.qvalue_cutoff)

