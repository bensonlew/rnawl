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
from mbio.api.database.ref_rna_v2.api_base import ApiBase
from mbio.packages.rna.dendrogram2newick import convert


class Wgcna(ApiBase):
    def __init__(self, bind_object):
        super(Wgcna, self).__init__(bind_object)
        self._project_type = 'tool_lab'

    def add_prepare_detail(self, prepare_result_dir, main_id):
        detail = dict()
        with open(prepare_result_dir+'/ignored_gene.list') as f:
            after_filter, before_filter = f.readline().strip('#').strip().split('/')
            detail.update(dict(after_filter=int(after_filter), before_filter=int(before_filter)))
        sample_cluster = glob.glob(prepare_result_dir + '/sample.cluster.dendrogram.txt')[0]
        linkage_pd = pd.read_table(sample_cluster, header=None)
        linkage_list = [list(x) for x in linkage_pd.as_matrix()]
        ordered_leaf = pd.read_table(sample_cluster[:-3] + 'order.txt', header=None)[0]
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
        independence = map(lambda x, y: x*y, list(scale_free['SFT.R.sq']), sploe)
        connectivity = list(scale_free['mean.k.'])
        detail.update(dict(power_list=powers, independence_list=independence, connectivity_list=connectivity))
        # create detail table
        detail.update(dict(prepare_id=ObjectId(main_id)))
        self.create_db_table('sg_wgcna_prepare_detail', [detail])
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
        # module size detail
        size_stat = glob.glob(result_dir + '/module_size.stat.xls')[0]
        size_stat_pd = pd.read_table(size_stat, header=0)
        size_stat_pd['pipeline_id'] = ObjectId(main_id)
        size_stat_pd['category'] = size_stat_pd['module']
        colour_list = size_stat_pd['module'].tolist()
        size_stat_pd['type'] = 'column'
        size_stat_dict = size_stat_pd.to_dict("records")
        self.create_db_table('sg_wgcna_module_stat_detail', size_stat_dict)
        dict_a = {"name": "module", "data": "size", "category": "category", "condition": {"type": "column"}}
        dict_b = json.dumps(dict_a, sort_keys=True, separators=(',', ':'))
        self.update_db_record('sg_wgcna_pipeline', ObjectId(main_id), column_data_detail_module=dict_b, colour_list_module=colour_list)

        # membership detail
        membership = glob.glob(result_dir + '/membership.xls')[0]
        membership_pd = pd.read_table(membership, header=0)
        # membership_pd.columns = ['seq_id', 'module', 'kme', 'block_id']
        membership_pd.set_index('seq_id', inplace=True)
        exp_pd = pd.read_table(exp_matrix, header=0, index_col=0)
        exp_pd.index.name = 'seq_id'
        membership_pd = membership_pd.join(g2n)
        membership_pd.to_csv(result_dir+"/gene_module_detail.xls", header=True, index=True, sep='\t')

        # sg_wgcna_module_cluster_detail
        with open(result_dir + '/module_corr.tree.txt') as f:
            tree_str = f.read().strip()
            labels = re.findall('[(,]([^(]*?):', tree_str)
            module_cluster_data_h = dict(
                data=tree_str,
                labels=labels,
                pipeline_id=ObjectId(main_id),
                type='tree',
                direction='h',
                group=[{'groupname': i, 'value': [i]} for i in labels]
            )
            module_cluster_data_v = dict(
                data=tree_str,
                labels=labels,
                pipeline_id=ObjectId(main_id),
                type='tree',
                direction='v',
                group=[{'groupname': i, 'value': [i]} for i in labels]
            )
        self.create_db_table('sg_wgcna_module_cluster_detail', [module_cluster_data_h])
        self.create_db_table('sg_wgcna_module_cluster_detail', [module_cluster_data_v])
        dict_a = {"name": "name", "condition": {"type": "tree"}}
        dict_b = json.dumps(dict_a, sort_keys=True, separators=(',', ':'))
        self.update_db_record('sg_wgcna_pipeline', ObjectId(main_id), cluster_data_detail_module = dict_b)
        heatmap_bar_list_h = list()
        for j in labels:
            heatmap_bar = dict(
                name=j,
                direction='h',
                type='heatmap_bar',
                level_color=j,
                pipeline_id=ObjectId(main_id)
            )
            heatmap_bar_list_h.append(heatmap_bar)
        heatmap_bar_list_v = list()
        for j in labels:
            heatmap_bar = dict(
                name=j,
                direction='v',
                type='heatmap_bar',
                level_color=j,
                pipeline_id=ObjectId(main_id)
            )
            heatmap_bar_list_v.append(heatmap_bar)
        heatmap_bar_data_h = dict(name = 'name', condition = {'type': 'heatmap_bar'})
        heatmap_bar_data_h = json.dumps(heatmap_bar_data_h)
        heatmap_bar_data_v = dict(name = 'name', condition = {'type': 'heatmap_bar'})
        heatmap_bar_data_v = json.dumps(heatmap_bar_data_v)
        self.create_db_table('sg_wgcna_heatmap_bar', heatmap_bar_list_h)
        self.create_db_table('sg_wgcna_heatmap_bar', heatmap_bar_list_v)
        self.update_db_record('sg_wgcna_pipeline', ObjectId(main_id), heatmap_bar_data_h=heatmap_bar_data_h)
        self.update_db_record('sg_wgcna_pipeline', ObjectId(main_id), heatmap_bar_data_v=heatmap_bar_data_v)
        self.update_db_record('sg_wgcna_pipeline', ObjectId(main_id), labels=labels)

        # sg_wgcna_module_corr_detail
        module_corr = glob.glob(result_dir + '/module_corr.matrix.xls')[0]
        module_corr_pd = pd.read_table(module_corr, header=0, index_col=0)
        module_corr_pd.index.name = "module"
        module_corr_pd = module_corr_pd.loc[labels, :]
        module_corr_pd.reset_index(inplace=True)
        module_corr_dict_list = module_corr_pd.to_dict('records')
        module_corr_dict_list = self.order_row_dict_list(module_corr_dict_list, labels + ['module'])
        self.create_db_table('sg_wgcna_module_corr_detail', module_corr_dict_list, tag_dict=dict(pipeline_id=ObjectId(main_id), type='heatmap'))
        dict_a = {"name": "module", "data": labels, "condition": {"type": "heatmap"}}
        dict_b = json.dumps(dict_a, sort_keys=True, separators=(',', ':'))
        self.update_db_record('sg_wgcna_pipeline', ObjectId(main_id), heatmap_data_detail_module = dict_b)
        # # sg_wgcna_module_eigengenes_detail
        # eigengenes = glob.glob(result_dir + '/eigengenes.txt')[0]
        # eigengenes_pd = pd.read_table(eigengenes, header=0, index_col=0)
        # module_list = list(eigengenes_pd.index)
        # sample_list = list(eigengenes_pd.columns)
        # eigengenes_pd.index.name = 'module'
        # eigengenes_pd.reset_index(inplace=True)
        # eigengenes_dict_list = eigengenes_pd.to_dict("records")
        # eigengenes_dict_list = self.order_row_dict_list(eigengenes_dict_list, list(eigengenes_pd.columns) + ['module'])
        # self.create_db_table('sg_wgcna_module_eigengenes_detail', eigengenes_dict_list, tag_dict=dict(module_id=ObjectId(main_id)))

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
        final_corr['pipeline_id'] = ObjectId(main_id)
        corr_all_dict_list = final_corr.to_dict("records")
        self.create_db_table('sg_wgcna_relate_gene_detail',corr_all_dict_list)

        # get gene significant stat
        # corr_all = corr_pd.abs().join(seq_annot_pd.loc[:, ["module"]])
        # mean_gs = corr_all.groupby("module").mean()
        # gs_std = corr_all.groupby("module").std()
        # gs_all = mean_gs.join(gs_std, lsuffix="_gs", rsuffix="_gs_std")
        # gs_all.index = ["ME"+x for x in gs_all.index]

        # module_trait.correlation.xls, module_trait.correlation_pvalues.xls
        module_trait_corr = glob.glob(result_dir + '/module_trait.correlation.xls')[0]
        module_trait_corr_pvalue = glob.glob(result_dir + '/module_trait.correlation_pvalues.xls')[0]
        corr = pd.read_table(module_trait_corr, index_col=0, header=0)
        # module_list = list(corr.index)
        trait_list = list(corr.columns)
        trait_list.insert(0, 'name')
        # trait_list_pvalue = [i + '_pvalue' for i in list(corr.columns)]
        corr_p = pd.read_table(module_trait_corr_pvalue, index_col=0, header=0)
        # corr_all = corr.join(corr_p, lsuffix='_corr', rsuffix="_pvalue")
        # corr_all.index.name = "module"
        corr.index.name = "module"
        corr_p.index.name = "module"
        module_size_pd = seq_annot_pd.groupby("module").count().iloc[:, [0]]
        module_size_pd.columns = ["size"]
        module_size_pd.index = ["ME"+x for x in module_size_pd.index]
        # corr_all = corr_all.join(module_size_pd)
        corr = corr.join(module_size_pd)
        corr_p = corr_p.join(module_size_pd)
        # corr_all = corr_all.join(gs_all)
        # corr_all['pipeline_id'] = ObjectId(main_id)
        corr['pipeline_id'] = ObjectId(main_id)
        corr_p['pipeline_id'] = ObjectId(main_id)
        # corr_all['type'] = 'heatmap'
        corr['type'] = 'heatmap'
        corr_p['type'] = 'heatmap_asterisk'
        colour_list = [i.lstrip('ME') for i in corr.index.tolist()]
        # corr_all.reset_index(inplace=True)
        corr.reset_index(inplace=True)
        corr_p.reset_index(inplace=True)
        # corr_all_dict_list = corr_all.to_dict("records")
        corr_dict_list = corr.to_dict("records")
        corr_p_dict_list = corr_p.to_dict("records")
        # self.create_db_table('sg_wgcna_relate_module_detail', corr_all_dict_list)
        self.create_db_table('sg_wgcna_relate_module_detail', corr_dict_list)
        self.create_db_table('sg_wgcna_relate_module_detail', corr_p_dict_list)
        relate_bar = dict(zip(corr['module'].tolist(), corr['size'].tolist()))
        relate_tree = dict(
            data="",
            pipeline_id=ObjectId(main_id),
            type='tree',
            direction='v',
            group=[{'groupname': i, 'value': [int(relate_bar[i])]} for i in relate_bar]
        )
        self.create_db_table('sg_wgcna_relate_module_tree', [relate_tree])

        relate_heatmap_bar_list_h = list()
        for j in corr['module'].tolist():
            heatmap_bar = dict(
                name=relate_bar[j],
                direction='v',
                type='heatmap_bar',
                level_color=j,
                pipeline_id=ObjectId(main_id)
            )
            relate_heatmap_bar_list_h.append(heatmap_bar)
        self.create_db_table('sg_wgcna_relate_module_bar', relate_heatmap_bar_list_h)
        dict_a = {"name": "module", "data": trait_list, "condition": {"type": "heatmap"}}
        dict_b = json.dumps(dict_a, sort_keys=True, separators=(',', ':'))
        dict_c = {'name': 'module', "data": trait_list, "condition": {"type": "heatmap_asterisk"}}
        dict_d = json.dumps(dict_c, sort_keys=True, separators=(',', ':'))
        dict_e = {"name": "name", "condition": {"type": "tree"}}
        dict_f = json.dumps(dict_e, sort_keys=True, separators=(',', ':'))
        relate_heatmap_bar_data_h = dict(name = 'name', condition = {'type': 'heatmap_bar'})
        relate_heatmap_bar_data_h = json.dumps(relate_heatmap_bar_data_h)
        self.update_db_record('sg_wgcna_pipeline', ObjectId(main_id), heatmap_data_detail_relate=dict_b,
                              colour_list_relate=colour_list, heatmap_asterisk_data=dict_d, heatmap_reladte_tree=dict_f,
                              heatmap_reladte_bar=relate_heatmap_bar_data_h
                              )


    def add_network_detail(self, result_dir, main_id, module):
        # for node
        if not glob.glob(result_dir+"/*network.nodes.txt"):
            self.update_db_record('sg_wgcna_network', main_id, status="end")
            print("!!! Network is Empty!")
            return
        node_file = glob.glob(result_dir+"/*network.nodes.txt")[0]
        node_pd = pd.read_table(node_file, header=0)
        node_pd['pipeline_id'] = ObjectId(main_id)
        node_pd['type'] = 'node'
        node_pd['module'] = module
        node_list_index = node_pd['id'].tolist()
        node_pd['node_id'] = range(0, len(node_list_index))
        node_dict_list = node_pd.to_dict("records")
        self.create_db_table("sg_wgcna_network_node_detail", node_dict_list)
        dict_a = {"name": "name", 'id': 'id', 'group': 'module', "condition": {"type": "node"}}
        dict_b = json.dumps(dict_a, sort_keys=True, separators=(',', ':'))
        self.update_db_record('sg_wgcna_pipeline', ObjectId(main_id), node_data = dict_b)
        # for edge
        edge_file = glob.glob(result_dir+"/*network.edges.txt")[0]
        edge_pd = pd.read_table(edge_file, header=0)
        edge_pd['pipeline_id'] = ObjectId(main_id)
        edge_pd['module'] = module
        edge_pd['source'] = [node_list_index.index(i) for i in edge_pd['source']]
        edge_pd['target'] = [node_list_index.index(i) for i in edge_pd['target']]
        edge_pd['name'] = edge_pd['source'].map(str) + '_' + edge_pd['target'].map(str)
        edge_pd['type'] = 'link'
        edge_dict_list = edge_pd.to_dict("records")
        self.create_db_table("sg_wgcna_network_edge_detail", edge_dict_list)
        dict_a = {"name": "name", "value": 'weight', "condition": {"type": "link"}}
        dict_b = json.dumps(dict_a, sort_keys=True, separators=(',', ':'))
        self.update_db_record('sg_wgcna_pipeline', ObjectId(main_id), link_data = dict_b)


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

