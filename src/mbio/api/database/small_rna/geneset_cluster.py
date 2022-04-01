#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""
@time    : 2018/11/7 17:51
@file    : geneset_cluster.py
"""

from mbio.api.database.small_rna.api_base import ApiBase
import json
import pandas as pd
import datetime
from bson.objectid import ObjectId
import os


class GenesetCluster(ApiBase):
    def __init__(self, bind_object):
        super(GenesetCluster, self).__init__(bind_object)

    def add_geneset_cluster(self, cluster_output_dir, main_id=None, project_sn='small_rna', task_id='dev_small_rna',
                            params=None):
        # prepare main_table data
        # print "type", type(main_id)
        # print "main_id", main_id
        results = os.listdir(cluster_output_dir)
        gene_cluster, sample_cluster = False, False
        genes, samples = list(), list()
        gene_tree, sample_tree = "", ""
        if "seq.cluster_tree.txt" in results:
            target_file = os.path.join(cluster_output_dir, "seq.cluster_tree.txt")
            with open(target_file) as f:
                gene_cluster = True
                gene_tree = f.readline().strip()
                genes = f.readline().strip().split(";")
        #
        if "sample.cluster_tree.txt" in results:
            target_file = os.path.join(cluster_output_dir, "sample.cluster_tree.txt")
            with open(target_file) as f:
                sample_cluster = True
                sample_tree = f.readline().strip()
                samples = f.readline().strip().split(";")
        #
        if "seq.kmeans_cluster.txt" in results:
            gene_cluster = True
            target_file = os.path.join(cluster_output_dir, "seq.kmeans_cluster.txt")
            with open(target_file) as f:
                genes = list()
                for line in f:
                    if not line.strip():
                        continue
                    genes += line.strip().split('\t')[1].split(";")
        #
        detail_info = list()
        trend_dict = dict()
        if ("seq.cluster_tree.txt" in results) or ("seq.kmeans_cluster.txt" in results) :
            sub_clusters = [x for x in results if x.startswith('seq.subcluster')]
            number_order = [(x, int(x.split('_')[1])) for x in sub_clusters]
            tmp = sorted(number_order, key=lambda x: x[1])
            sub_clusters = [x[0] for x in tmp]
            for sub in sub_clusters:
                target_file = os.path.join(cluster_output_dir, sub)
                tmp_df = pd.read_table(target_file, header=0)
                sub_cluster_id = int(sub.split('_')[1])
                tmp_df["sub_cluster"] = sub_cluster_id
                detail_info += json.loads(tmp_df.to_json(orient="records"))
                mean_dict = tmp_df.iloc[:, 1:-1].mean().to_dict()
                trend_dict[str(sub_cluster_id)] = mean_dict
        #
        target_file = os.path.join(cluster_output_dir, "expression_matrix.xls")
        exp_pd = pd.read_table(target_file, header=0)
        if not detail_info:
            detail_info = exp_pd.to_dict('records')
        if not genes:
            if 'seq_id' in exp_pd.columns:
                genes = list(exp_pd['seq_id'])
            elif 'accession_id' in exp_pd.columns:
                genes = list(exp_pd['accession_id'])
        if not samples:
            samples = list(exp_pd.columns)[1:]
        # add main table info'
        if main_id is None:
            name = "GeneSet_Cluster" + '_'
            time_now = datetime.datetime.now()
            name += time_now.strftime("%Y%m%d_%H%M%S")
            if params is None:
                params_dict = dict()
            elif type(params) == dict:
                params_dict = params
                params = json.dumps(params, sort_keys=True, separators=(',', ':'))
            else:
                params_dict = json.loads(params)

            main_info = dict(
                project_sn=project_sn,
                task_id=task_id,
                name=name,
                created_ts=time_now.strftime('%Y-%m-%d %H:%M:%S'),
                desc='geneset cluster main table',
                status="start",
                params=params,
            )
            main_id = self.create_db_table('sg_geneset_cluster', [main_info])
        else:
            try:
                main_id = ObjectId(main_id)
            except Exception as e:
                self.bind_object.set_error('main_id 格式错误')
        # update main table
        self.update_db_record('sg_geneset_cluster', main_id,
                              trend_dict=trend_dict,
                              samples=samples,
                              gene_cluster=gene_cluster,
                              sample_cluster=sample_cluster, )
        # add detail info
        tree_info = dict(
            genes=genes,
            gene_tree=gene_tree,
            sample_tree=sample_tree,
            cluster_id=main_id,
        )

        def tran_dic_key(dic_con):
            if not isinstance(dic_con, dict):
                dic_con = dict(dic_con)
            dic_con['seq_id'] = dic_con.pop('accession_id')
            return dic_con

        detail_info = [tran_dic_key(dic) for dic in detail_info]
        self.create_db_table('sg_geneset_cluster_tree', [tree_info])
        self.create_db_table('sg_geneset_cluster_detail', detail_info, tag_dict=dict(cluster_id=main_id))
        self.update_db_record('sg_geneset_cluster', main_id, status="end", main_id=main_id, )
        return main_id

    # def add_geneset_cluster(self, cluster_output_dir, gene_detail_file, main_id=None, project_sn='small_rna',
    #                         task_id='dev_small_rna',
    #                         params=None):
    #     # prepare main_table data
    #     results = os.listdir(cluster_output_dir)
    #     gene_cluster, sample_cluster = False, False
    #     genes, samples = list(), list()
    #     gene_tree, sample_tree = "", ""
    #
    #     if "seq.cluster_tree.txt" in results:
    #         target_file = os.path.join(cluster_output_dir, "seq.cluster_tree.txt")
    #         with open(target_file) as f:
    #             gene_cluster = True
    #             gene_tree = f.readline().strip()
    #             genes = f.readline().strip().split(";")
    #     #
    #     if "sample.cluster_tree.txt" in results:
    #         target_file = os.path.join(cluster_output_dir, "sample.cluster_tree.txt")
    #         with open(target_file) as f:
    #             sample_cluster = True
    #             sample_tree = f.readline().strip()
    #             samples = f.readline().strip().split(";")
    #     #
    #     if "seq.kmeans_cluster.txt" in results:
    #         gene_cluster = True
    #         target_file = os.path.join(cluster_output_dir, "seq.kmeans_cluster.txt")
    #         with open(target_file) as f:
    #             genes = list()
    #             for line in f:
    #                 if not line.strip():
    #                     continue
    #                 genes += line.strip().split('\t')[1].split(";")
    #     #
    #     detail_info = list()
    #     trend_dict = dict()
    #     seq_id2name = dict()
    #     with open(gene_detail_file, 'rb') as f:
    #         for linen in f.readlines()[1:]:
    #             line = linen.strip()
    #             if len(line.split(b"\t")) > 2:
    #                 if line.split(b"\t")[2].strip() in ["-", "_"]:
    #                     seq_id2name[line.split(b"\t")[0]] = line.split(b"\t")[0]
    #                     seq_id2name[line.split(b"\t")[1]] = line.split(b"\t")[1]
    #                 else:
    #                     seq_id2name[line.split(b"\t")[0]] = line.split(b"\t")[2].strip()
    #                     seq_id2name[line.split(b"\t")[1]] = line.split(b"\t")[2].strip()
    #             else:
    #                 seq_id2name[line.split(b"\t")[0]] = line.split(b"\t")[0]
    #                 seq_id2name[line.split(b"\t")[1]] = line.split(b"\t")[1]
    #
    #     if ("seq.cluster_tree.txt" in results) or ("seq.kmeans_cluster.txt" in results):
    #         sub_clusters = [x for x in results if x.startswith('seq.subcluster')]
    #         number_order = [(x, int(x.split('_')[1])) for x in sub_clusters]
    #         tmp = sorted(number_order, key=lambda x: x[1])
    #         sub_clusters = [x[0] for x in tmp]
    #         for sub in sub_clusters:
    #             target_file = os.path.join(cluster_output_dir, sub)
    #             tmp_df = pd.read_table(target_file, header=0)
    #             sub_cluster_id = int(sub.split('_')[1])
    #             tmp_df["sub_cluster"] = sub_cluster_id
    #             tmp_df["gene_name"] = tmp_df['seq_id'].map(lambda x: seq_id2name[x])
    #             detail_info += json.loads(tmp_df.to_json(orient="records"))
    #             mean_dict = tmp_df.iloc[:, 1:-1].mean().to_dict()
    #             trend_dict[str(sub_cluster_id)] = mean_dict
    #     #
    #     target_file = os.path.join(cluster_output_dir, "expression_matrix.xls")
    #
    #     exp_pd = pd.read_table(target_file, header=0)
    #
    #     if not detail_info:
    #         exp_pd['gene_name'] = exp_pd['seq_id'].map(lambda x: seq_id2name[x])
    #         detail_info = exp_pd.to_dict('records')
    #     if not genes:
    #         genes = list(exp_pd['seq_id'])
    #     if not samples:
    #         samples = list(exp_pd.columns)[1:]
    #     # add main table info'
    #     if main_id is None:
    #         name = "GeneSet_Cluster" + '_'
    #         time_now = datetime.datetime.now()
    #         name += time_now.strftime("%Y%m%d_%H%M%S")
    #         if params is None:
    #             params_dict = dict()
    #         elif type(params) == dict:
    #             params_dict = params
    #             params = json.dumps(params, sort_keys=True, separators=(',', ':'))
    #         else:
    #             params_dict = json.loads(params)
    #
    #         main_info = dict(
    #             project_sn=project_sn,
    #             task_id=task_id,
    #             name=name,
    #             created_ts=time_now.strftime('%Y-%m-%d %H:%M:%S'),
    #             desc='geneset cluster main table',
    #             status="start",
    #             params=params,
    #             type='T' if "exp_level" not in params_dict else params_dict["exp_level"],
    #         )
    #         main_id = self.create_db_table('sg_geneset_cluster', [main_info])
    #     else:
    #         if isinstance(main_id, types.StringTypes):
    #             main_id = ObjectId(main_id)
    #     # update main table
    #     self.update_db_record('sg_geneset_cluster', main_id,
    #                           trend_dict=trend_dict,
    #                           samples=samples,
    #                           gene_cluster=gene_cluster,
    #                           sample_cluster=sample_cluster, )
    #     # add detail info
    #     tree_info = dict(
    #         genes=genes,
    #         gene_tree=gene_tree,
    #         sample_tree=sample_tree,
    #         cluster_id=main_id,
    #     )
    #     self.create_db_table('sg_geneset_cluster_tree', [tree_info])
    #     self.create_db_table('sg_geneset_cluster_detail', detail_info, tag_dict=dict(cluster_id=main_id))
    #     self.update_db_record('sg_geneset_cluster', main_id, status="end", main_id=main_id, )
    #     return main_id
