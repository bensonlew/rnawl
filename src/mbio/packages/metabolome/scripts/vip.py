# -*- coding: utf-8 -*-
# __author__ = 'shaohua.yuan'
# last_modify:20180615

import fastcluster as hclust
import sys
import argparse
import os, copy
import pandas as pd
import shutil
import scipy
from scipy import stats
from sklearn import decomposition, preprocessing
from mbio.packages.metabolome.scripts.corr_cluster import CorrClust
import numpy as np

def vip_filter(diff_profile, outDir, vip_type, vip_cut, vip_top, s_cluster, g_cluster, c_clu_type="hierarchy",
               r_clu_type="hierarchy", c_clu_method="complete", r_clu_method="complete", c_dist_method="euclidean",
               r_dist_method="euclidean", n_clusters=0):
    """
    Vip分析，根据阈值筛选top vip
    """
    diff_table = pd.read_table(diff_profile, sep="\t", header=0, index_col=0)
    if vip_type == "plsda":
        vip_col = "Vip_plsda"
    else:
        vip_col = "Vip_oplsda"
    filter = diff_table.loc[diff_table[vip_col] > float(vip_cut)]
    top = filter.sort_values([vip_col], ascending=False)[0:int(vip_top)]
    outfile = outDir + "/Vip_tmp.xls"
    #outfile = outDir + "/Vip.xls"
    top.to_csv(outfile, sep="\t")
    cluster = CorrClust()
    if len(top) < 2:
        g_cluster = None
    cluster.cor_cluster("clu", top, s_cluster, g_cluster, c_clu_type, r_clu_type, c_clu_method, r_clu_method,
                        c_dist_method, r_dist_method, outDir, n_clusters=n_clusters)
    return top

def select_metab(exp_file, outDir, metab_select=None, group_file=None, group_name=None, group_method=None):
    exp_table = pd.read_table(exp_file, sep="\t", header=0)
    exp_table.rename(columns={exp_table.columns[0]: "metab_id"}, inplace=True)
    select_table = exp_table
    if metab_select:
        metab_list = pd.read_table(metab_select, sep='\t', header=None)
        metab_list.columns = ["metab_id"]
        print metab_list.head()
        select_table = pd.merge(exp_table, metab_list, how="inner", on="metab_id")
        print "metablsit"
        #print select_table.head()
    if group_file:
        samples = get_samples(group_file, group_name=group_name)
        select_columns = ["metab_id", "Vip_plsda", "Vip_oplsda", "P_value", "FC"] + samples
        select_table = select_table[select_columns]
        print "group"
        #print select_table.head()
    if group_method:
        select_columns = ["metab_id", "Vip_plsda", "Vip_oplsda", "P_value", "FC"]
        other = select_table[select_columns]
        group_dict = get_group_dict(group_file, group_name=group_name)
        if group_method == "sum":
            groups = select_table[group_dict.keys()].groupby(group_dict, axis=1).sum()
        elif group_method == "average":
            groups = select_table[group_dict.keys()].groupby(group_dict, axis=1).mean()
        elif group_method == "median":
            groups = select_table[group_dict.keys()].groupby(group_dict, axis=1).median()
        groups.index = select_table.index
        select_table = pd.merge(other, groups, how="inner", left_index=True, right_index=True)
    outfile = outDir + "/select_metab.xls"
    select_table.to_csv(outfile, sep="\t", index=False)
    return outfile

def get_samples(group_file, group_name=None):
    group_table = pd.read_table(group_file, sep="\t", header=0)
    if group_name:
        group_names = group_name.split(",")
        group_table = group_table[group_table.iloc[:, 1].isin(group_names)]
    sample_list = group_table.iloc[:, 0].tolist()
    return sample_list

def get_group_dict(group_file, group_name=None):
    group_dict = {}
    print group_name
    if group_name:
        group_names = group_name.split(",")
    with open(group_file) as f:
        f.next()
        for line in f:
            line = line.strip().split("\t")
            sample = line[0]
            group = line[1]
            if group_name:
                if group in group_names:
                    group_dict[sample] = group
            else:
                group_dict[sample] = group
    return group_dict

def scale_data(ori_table, all_abu, outDir, not_all=True):

    table = pd.read_table(ori_table, sep="\t", index_col=0)
    all_metab_list = ["Vip_plsda", "Vip_oplsda", "P_value", "FC", "fdr"]  ##new add fdr 20190618
    names = table.columns.tolist()
    metab_list= []
    for x in all_metab_list:
        if x in names:
            metab_list.append(x)
    metab_other = table.loc[:, metab_list]
    samples = copy.deepcopy(names)
    [samples.remove(i) for i in metab_list]

    all_samples_abu = pd.read_table(all_abu, sep="\t", index_col=0)
    if not_all == False:
        scaled_data = all_samples_abu.apply(lambda x: (x - np.mean(x)) / np.std(x, ddof=1), axis=1)
    else:   #标准化不用所有样本的数据  20190824 zouguanqing
        print('scale sample: %s' % ','.join(samples))
        sub_abu = all_samples_abu[samples]
        scaled_data = sub_abu.apply(lambda x: (x - np.mean(x)) / np.std(x, ddof=1), axis=1)

    sample_table = scaled_data.loc[:, samples]
    scaled_data = pd.merge(metab_other, sample_table, how="inner", left_index=True, right_index=True)
    exp_profile = os.path.join(outDir, "scale_data.xls")
    scaled_data.to_csv(exp_profile, index=True, header=True, sep="\t")
    return exp_profile

def orgin_abu_from_scale(origin_file, scale_top_file, outDir, group_method=None, group_file=None, group_name=None):
    table = pd.read_table(scale_top_file, sep="\t", index_col=0)
    select_names = table.index
    origin_table = pd.read_table(origin_file, sep="\t", index_col=0)
    selecl_origin = origin_table.loc[select_names,]
    select_abu = selecl_origin
    all_metab_list = ["Vip_plsda", "Vip_oplsda", "P_value", "FC", "fdr"] ##new add fdr 20190618
    names_col = select_abu.columns.tolist()
    metab_list= []
    for x in all_metab_list:
        if x in names_col:
            metab_list.append(x)

    metab_other = select_abu.loc[:, metab_list]
    samples = copy.deepcopy(names_col)
    [samples.remove(i) for i in metab_list]
    sample_table = select_abu.loc[:, samples]
    if group_method:
        if not group_file:
            raise Exception("选择分组计算时没有分组文件")
        dict_group = get_group_dict(group_file, group_name=group_name)
        if group_method == "sum":
            groups = sample_table[dict_group.keys()].groupby(dict_group, axis=1).sum()
        elif group_method == "average":
            groups = sample_table[dict_group.keys()].groupby(dict_group, axis=1).mean()
        elif group_method == "median":
            groups = sample_table[dict_group.keys()].groupby(dict_group, axis=1).median()
        groups.index = sample_table.index
        select_abu = pd.merge(metab_other, groups, how="inner", left_index=True, right_index=True)
    select_outfile = os.path.join(outDir, "Vip_before_scale.xls")
    select_abu.to_csv(select_outfile, sep="\t", index=True, header=True)

if __name__ == '__main__':
    parser = argparse.ArgumentParser()
    parser.add_argument('-exp', type=str, metavar="exp_matrix_file", required=True,
                        help="expression matrix file")
    #parser.add_argument('-group', type=str, metavar="group_file", default=None, help="file with two col: sample\tgroup")
    parser.add_argument('-out', type=str, default=None, help="default is local dir. Output directory name")
    parser.add_argument('--T', action='store_true', help="need transform exp_file")
    parser.add_argument('-n_clusters', type=int, default=0, metavar="cluster_num", help="expected_cluster_number")
    parser.add_argument('--nsc', action='store_true', help="'no-sample-cluster', do not perform sample cluster")
    parser.add_argument('--ngc', action='store_true', help="'no-gene-cluster', do not perform gene cluster")
    parser.add_argument('-sct', metavar='sample-cluster-type', default='', help="hierarchy or kmeans cluster")
    parser.add_argument('-gct', metavar='gene-cluster-type', default='hierarchy', help="hierarchy or kmeans cluster")
    parser.add_argument('-scm', metavar='sample-cluster-method', default='complete', help="hclust method")
    parser.add_argument('-gcm', metavar='gene-cluster-method', default='complete', help="hclust method")
    parser.add_argument('-scd', metavar='sample-cluster-distance', default='correlation', help="hclust distance metric")
    parser.add_argument('-gcd', metavar='gene-cluster-distance', default='euclidean', help="hclust distance metric")
    parser.add_argument('-vty', metavar='vip_source', default='oplsda', help="vip from plsda or oplsda")
    parser.add_argument('-vcut', metavar='vip cutoff', default=1, help="vip cutoff")
    parser.add_argument('-vtop', metavar='top_vip', default=30, help="top_vip")
    parser.add_argument('-group', metavar='group file', help="group file")
    parser.add_argument('-dn', metavar='use group name', help="use group name in group file")
    parser.add_argument('-metaf', metavar='select metab file', help="select metab file")
    parser.add_argument('-gm', metavar='group method', default="none", help="group method")
    parser.add_argument('--scale', action='store_true', help="scale table")
    parser.add_argument('-t_abu', metavar='total_abu', help="all samples abu table for scale")
    args = parser.parse_args()
    exp_pd = args.exp
    outDir = args.out
    s_cluster, g_cluster = True, True
    trans = False
    scale = False
    metab_select, group_file, group_name, group_method = None, None, None, None
    myorder = CorrClust()
    myorder.remove_file(outDir)
    if args.nsc:
        s_cluster = False
    if args.ngc:
        g_cluster = False
    if args.metaf:
        metab_select = args.metaf
    if args.gm != "none":
        group_file = args.group
        group_method = args.gm
    if args.dn:
        diff_group_name = args.dn
    if args.scale:
        if not args.t_abu:
            raise Exception('scale时必须输入-t_abu所有样本的丰度表!')
        scale = True
        exp_scale = scale_data(exp_pd, args.t_abu, outDir)
        exp_pd = select_metab(exp_scale, outDir, metab_select=metab_select, group_file=group_file,
                              group_name=diff_group_name, group_method=group_method)
    else:
        exp_pd = select_metab(exp_pd, outDir, metab_select=metab_select, group_file=group_file,
                              group_name=diff_group_name, group_method=group_method)
    tmp_table = vip_filter(exp_pd, outDir, args.vty, args.vcut, args.vtop, s_cluster, g_cluster, c_clu_type=args.sct,
                           r_clu_type=args.gct, c_clu_method=args.scm, r_clu_method=args.gcm, c_dist_method=args.scd,
                           r_dist_method=args.gcd, n_clusters=args.n_clusters)
    col_tree_file = None
    row_tree_file = None
    if not args.nsc and args.sct == "hierarchy":
        col_tree_file = outDir + "/col.cluster_tree.txt"
    if not args.ngc and args.gct == "hierarchy" and len(tmp_table) > 1:
        row_tree_file = outDir + "/row.cluster_tree.txt"
    if args.gct == "hierarchy":
        myorder.order_sample(tmp_table, outDir + "/Vip.xls", coltree=col_tree_file, rowtree=row_tree_file)
    elif args.gct == "kmeans":
        if len(tmp_table) > 1:
            myorder.order_kmeans(tmp_table, outDir + "/Vip.xls", outDir, cluster=True)
        else:
            myorder.order_kmeans(tmp_table, outDir + "/Vip.xls", outDir)
    else:
        os.link(tmp_file, outDir + "/Vip.xls")
    if args.scale:
        order_file = outDir + "/Vip.xls"
        orgin_abu_from_scale(args.exp, order_file, outDir, group_method, group_file=group_file,
                             group_name=diff_group_name)
