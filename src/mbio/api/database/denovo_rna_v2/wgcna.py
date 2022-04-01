# !/usr/bin/python
# -*- coding: utf-8 -*-
from bson.objectid import ObjectId
import types
import re
import os
import json
import pandas as pd
import numpy as np
from scipy import stats
import math
from sklearn import decomposition
import scipy.cluster.hierarchy as sch
import fastcluster as hclust
from collections import OrderedDict
import unittest
import datetime
import json
import glob
from api_base import ApiBase


class Wgcna(ApiBase):
    def __init__(self, bind_object):
        super(Wgcna, self).__init__(bind_object)

    def add_prepare_detail(self, prepare_result_dir, main_id):
        detail = dict()
        with open(prepare_result_dir+'/ignored_gene.list') as f:
            after_filter, before_filter = f.readline().strip('#').strip().split('/')
            detail.update(dict(after_filter=int(after_filter), before_filter=int(before_filter)))
        # add sample cluster info
        # with open(prepare_result_dir+'/dendrogram.json') as f:
        #     dendrogram_dict = json.loads(f.readline().strip())
        #     min_h = float(f.readline().strip())
        #     detail.update(dendrogram_dict)
        sample_cluster = glob.glob(prepare_result_dir + '/sample.cluster.dendrogram.txt')[0]
        linkage_pd = pd.read_table(sample_cluster, header=None)
        linkage_list = [list(x) for x in linkage_pd.as_matrix()]
        ordered_leaf = pd.read_table(sample_cluster[:-3] + 'order.txt', header=None)[0]
        # x_coord, y_coord, _ = self.get_tree_coord(linkage_pd, ordered_leaf, height_ratio=0.06)
        labels = [x.strip() for x in open(sample_cluster[:-3] + 'labels.txt')]
        detail.update(dict(
            linkage=linkage_list,
            order=list(ordered_leaf),
            ivl=labels,
        ))
        # power select
        power_record = glob.glob(prepare_result_dir+'/powerEstimate_*')[0]
        power_estimate = os.path.basename(power_record).split('_')[1]
        if power_estimate == "NA":
            power_estimate = 6
        else:
            power_estimate = int(power_estimate)
        detail.update(dict(power_estimate=power_estimate))
        # scale_free analysis
        scale_free_analysis = prepare_result_dir+'/scale_free_analysis.xls'
        scale_free = pd.read_table(scale_free_analysis, index_col=0, header=0)
        powers = list(scale_free.index)
        sploe = [-1 if i > 0 else 1 for i in list(scale_free['slope'])]
        independence = map(lambda x, y: x * y, list(scale_free['SFT.R.sq']), sploe)
        connectivity = list(scale_free['mean.k.'])
        detail.update(dict(power_list=powers, independence_list=independence, connectivity_list=connectivity))
        # create detail table
        detail.update(dict(prepare_id=ObjectId(main_id)))
        self.create_db_table('sg_wgcna_prepare_detail', [detail])
        self.update_db_record('sg_wgcna_prepare', main_id, status="end", )
        return power_estimate

    @staticmethod
    def get_tree_coord(linkage_pd, ordered_leaf, height_ratio=0.01):
        """
        format block tree data
        :param ordered_leaf: [1,2,3,..],代表画聚类树的结果
        :param linkage_pd: [n1,n2,height] 来自R的聚类结果，不适用python的聚类结果
        :return:
        """
        sn = len(ordered_leaf)
        number = sn
        leaf_height_dict = {x: 0 for x in range(sn)}
        origin_leaf_x_loc = range(5,(sn - 1) * 10 + 6,10)
        leaf_x_dict = dict(zip(ordered_leaf, origin_leaf_x_loc))
        x_coord = list()
        y_coord = list()
        steps = list()
        mean_height = linkage_pd.iloc[:,2].mean()
        for each in linkage_pd.iterrows():
            sn += 1
            n1, n2 = int(each[1][0]), int(each[1][1])
            if n1 > 0:
                n1 = number + n1
            else:
                n1 = abs(n1)
            if n2 > 0:
                n2 = number + n2
            else:
                n2 = abs(n2)
            dis = float(each[1][2])
            leaf_height_dict[sn] = dis
            leaf_x_dict[sn] = (leaf_x_dict[n1] + leaf_x_dict[n2]) / 2
            if n1 <= number:
                y1 = dis - mean_height * height_ratio
            else:
                y1 = leaf_height_dict[n1]
            y2, y3 = dis, dis
            if n2 <= number:
                y4 = dis - mean_height * height_ratio
            else:
                y4 = leaf_height_dict[n2]
            #     y1 = leaf_height_dict[n1]
            #     y4 = leaf_height_dict[n2]
            steps.append([n1, n2])
            y_coord.append([y1, y2, y3, y4])
            x_coord.append([leaf_x_dict[n1], leaf_x_dict[n1], leaf_x_dict[n2], leaf_x_dict[n2]])
        return x_coord,y_coord,ordered_leaf

    def add_module_detail(self, result_dir, gene_id2gene_name, exp_matrix, main_id):
        """
        模块识别结果入库
        :param result_dir:
        :param main_id:
        :return:
        """
        # gene_id2gene_name
        g2n = pd.read_table(gene_id2gene_name, header=0, index_col=0)
        g2n.index.name = "seq_id"

        # block dendrogram -> sg_wgcna_module_tree_detail
        blocks = glob.glob(result_dir + '/block_*.dendrogram.txt')
        blocks_data = list()
        for each in blocks:
            linkage_pd = pd.read_table(each, header=None)
            linkage_list = [list(x) for x in linkage_pd.as_matrix()]
            ordered_leaf = pd.read_table(each[:-3] + 'order.txt', header=None)[0]
            # x_coord, y_coord, _ = self.get_tree_coord(linkage_pd, ordered_leaf)
            labels = [x.strip() for x in open(each[:-3] + 'labels.txt')]
            ordered_seqs = [x.strip() for x in open(each[:-3] + 'ordered_seqs.txt')]
            blocks_data.append(dict(
                # icoord=x_coord,
                # dcoord=y_coord,
                linkage=linkage_list,
                ivl=labels,
                order=list(ordered_leaf),
                ordered_seqs=ordered_seqs,
                block_id=os.path.basename(each).split(".")[0].split("_")[1],
                module_id=ObjectId(main_id),
            ))
        self.create_db_table('sg_wgcna_module_tree_detail', blocks_data)

        # module size detail
        size_stat = glob.glob(result_dir + '/module_size.stat.xls')[0]
        size_stat_pd = pd.read_table(size_stat, header=0)
        size_stat_pd['module_id'] = ObjectId(main_id)
        size_stat_dict = size_stat_pd.to_dict("records")
        self.create_db_table('sg_wgcna_module_stat_detail', size_stat_dict)

        # membership detail
        membership = glob.glob(result_dir + '/membership.xls')[0]
        membership_pd = pd.read_table(membership, header=0)
        # membership_pd.columns = ['seq_id', 'module', 'kme', 'block_id']
        membership_pd.set_index('seq_id', inplace=True)
        exp_pd = pd.read_table(exp_matrix, header=0, index_col=0)
        exp_pd.index.name = 'seq_id'
        membership_pd = membership_pd.join(g2n)
        membership_pd.to_csv(result_dir+"/gene_module_detail.xls", header=True, index=True, sep='\t')
        membership_pd = membership_pd.join(exp_pd)
        membership_pd['module_id'] = ObjectId(main_id)
        membership_pd.reset_index(inplace=True)
        membership_dict_list = membership_pd.to_dict('records')
        self.create_db_table('sg_wgcna_module_membership_detail', membership_dict_list)

        # sg_wgcna_module_unmerged_detail
        module_cluster = glob.glob(result_dir + '/unmerged.module_corr.dendrogram.txt')[0]
        linkage_pd = pd.read_table(module_cluster, header=None)
        linkage_list = [list(x) for x in linkage_pd.as_matrix()]
        ordered_leaf = pd.read_table(module_cluster[:-3] + 'order.txt',header=None)[0]
        # x_coord, y_coord, _ = self.get_tree_coord(linkage_pd, ordered_leaf, height_ratio=0.06)
        labels = [x.strip() for x in open(module_cluster[:-3] + 'labels.txt')]
        unmerged_module_cluster_data = dict(
            # icoord=x_coord,
            # dcoord=y_coord,
            linkage=linkage_list,
            order=list(ordered_leaf),
            ivl=labels,
            module_id=ObjectId(main_id),
        )
        self.create_db_table('sg_wgcna_module_unmerged_detail',[unmerged_module_cluster_data])

        # sg_wgcna_module_cluster_detail
        with open(result_dir + '/module_corr.tree.txt') as f:
            tree_str = f.read().strip()
            labels = re.findall('[(,]([^(]*?):', tree_str)
            module_cluster_data = dict(
                cluster_tree=tree_str,
                labels=labels,
                module_id=ObjectId(main_id),
            )
        self.create_db_table('sg_wgcna_module_cluster_detail', [module_cluster_data])

        # sg_wgcna_module_corr_detail
        module_corr = glob.glob(result_dir + '/module_corr.matrix.xls')[0]
        module_corr_pd = pd.read_table(module_corr, header=0, index_col=0)
        # origin_colnames = module_corr_pd.columns
        # origin_rownames = module_corr_pd.index
        # module_corr_pd.columns = [x[2:] for x in origin_colnames]
        # module_corr_pd.index = [x[2:] for x in origin_rownames]
        module_corr_pd.index.name = "module"
        module_corr_pd = module_corr_pd.loc[labels, :]
        module_corr_pd.reset_index(inplace=True)
        module_corr_dict_list = module_corr_pd.to_dict('records')
        module_corr_dict_list = self.order_row_dict_list(module_corr_dict_list, labels + ['module'])
        self.create_db_table('sg_wgcna_module_corr_detail', module_corr_dict_list, tag_dict=dict(module_id=ObjectId(main_id)))

        # sg_wgcna_module_eigengenes_detail
        eigengenes = glob.glob(result_dir + '/eigengenes.txt')[0]
        eigengenes_pd = pd.read_table(eigengenes, header=0, index_col=0)
        module_list = list(eigengenes_pd.index)
        sample_list = list(eigengenes_pd.columns)
        eigengenes_pd.index.name = 'module'
        eigengenes_pd.reset_index(inplace=True)
        eigengenes_dict_list = eigengenes_pd.to_dict("records")
        eigengenes_dict_list = self.order_row_dict_list(eigengenes_dict_list, list(eigengenes_pd.columns) + ['module'])
        self.create_db_table('sg_wgcna_module_eigengenes_detail', eigengenes_dict_list, tag_dict=dict(module_id=ObjectId(main_id)))

        # update main table
        self.update_db_record('sg_wgcna_module', main_id, status="end", modules=module_list, samples=sample_list, block_num=len(blocks))

    def add_relate_detail(self, result_dir, seq_annot, main_id):
        """
        导入表型关联分析结果
        :param result_dir:
        :param seq_annot:
        :param main_id:
        :return:
        """
        # seq annot
        seq_annot_pd = pd.read_table(seq_annot, header=0)
        seq_annot_pd.set_index('seq_id', inplace=True)

        # gene_trait.correlation.xls == gene significant detail
        gene_trait_corr = glob.glob(result_dir + '/gene_trait.correlation.xls')[0]
        corr_pd = pd.read_table(gene_trait_corr,index_col=0,header=0)
        corr_pd.index.name = 'seq_id'
        corr_all = corr_pd.join(seq_annot_pd)
        final_corr = corr_all.reset_index()
        final_corr['relate_id'] = ObjectId(main_id)
        corr_all_dict_list = final_corr.to_dict("records")
        self.create_db_table('sg_wgcna_relate_gene_detail',corr_all_dict_list)

        # get gene significant stat
        corr_all = corr_pd.abs().join(seq_annot_pd.loc[:, ["module"]])
        mean_gs = corr_all.groupby("module").mean()
        gs_std = corr_all.groupby("module").std()
        gs_all = mean_gs.join(gs_std, lsuffix="_gs", rsuffix="_gs_std")
        gs_all.index = ["ME"+x for x in gs_all.index]

        # module_trait.correlation.xls, module_trait.correlation_pvalues.xls
        module_trait_corr = glob.glob(result_dir + '/module_trait.correlation.xls')[0]
        module_trait_corr_pvalue = glob.glob(result_dir + '/module_trait.correlation_pvalues.xls')[0]
        corr = pd.read_table(module_trait_corr, index_col=0, header=0)
        module_list = list(corr.index)
        trait_list = list(corr.columns)
        corr_p = pd.read_table(module_trait_corr_pvalue, index_col=0, header=0)
        corr_all = corr.join(corr_p, lsuffix='_corr', rsuffix="_pvalue")
        corr_all.index.name = "module"
        module_size_pd = seq_annot_pd.groupby("module").count().iloc[:, [0]]
        module_size_pd.columns = ["size"]
        module_size_pd.index = ["ME"+x for x in module_size_pd.index]
        corr_all = corr_all.join(module_size_pd)
        corr_all = corr_all.join(gs_all)
        corr_all['relate_id'] = ObjectId(main_id)
        corr_all.reset_index(inplace=True)
        corr_all_dict_list = corr_all.to_dict("records")
        self.create_db_table('sg_wgcna_relate_module_detail', corr_all_dict_list)
        # update main table
        self.update_db_record('sg_wgcna_relate', main_id, status="end", modules=module_list, traits=trait_list)

    def add_network_detail(self, result_dir, main_id):
        # for node
        if not glob.glob(result_dir+"/*network.nodes.txt"):
            self.update_db_record('sg_wgcna_network', main_id, status="end")
            print("!!! Network is Empty!")
            return
        node_file = glob.glob(result_dir+"/*network.nodes.txt")[0]
        node_pd = pd.read_table(node_file, header=0)
        node_pd['network_id'] = ObjectId(main_id)
        node_dict_list = node_pd.to_dict("records")
        self.create_db_table("sg_wgcna_network_node_detail", node_dict_list)
        # for edge
        edge_file = glob.glob(result_dir+"/*network.edges.txt")[0]
        edge_pd = pd.read_table(edge_file, header=0)
        edge_pd['network_id'] = ObjectId(main_id)
        edge_dict_list = edge_pd.to_dict("records")
        self.create_db_table("sg_wgcna_network_edge_detail", edge_dict_list)
        # update main table
        self.update_db_record('sg_wgcna_network', main_id, status="end")

    def remove_all_test_result(self):
        target_coll = list()
        for each in self.db.collection_names():
            if each.startswith("sg_wgcna_"):
                target_coll.append(each)
        for each in target_coll:
            self.db.drop_collection(each)


if __name__ == '__main__':
    wgcna = Wgcna(None)
    wgcna.add_module_detail("/mnt/ilustre/users/sanger-dev/sg-users/deqing/workspace/20180319/WgcnaModule_tsg_28226_50578_372539/WgcnaModule"
                            ,"/mnt/ilustre/users/sanger-dev/workspace/20180319/WgcnaModule_tsg_28226_50578_372539/seq_id2gene_name.txt",
                            "5aaf7efca4e1af606ad5e527")

