# -*- coding: utf-8 -*-
# __author__ = 'zoujiaxun'
# last_modify:20200812
from biocluster.api.database.base import Base, report_check
from bson.objectid import ObjectId
from types import StringTypes
from bson.son import SON
from biocluster.config import Config
import datetime
import json
import os
import pandas as pd
import glob
import re
from mbio.api.database.ref_rna_v2.api_base import ApiBase

class Gsva(ApiBase):
    def __init__(self, bind_object):
        super(Gsva, self).__init__(bind_object)
        self._project_type = 'tool_lab'

    def add_gsva(self, main_id, es_matrix, gmx, params=None, project_sn='medical_transcriptome', task_id='medical_transcriptome'):
        """
        dump exp data into database
        :param exp_matrix: expression matrix path or an express matrix in pandas DataFrame format.
        :param exp_level: str, transcript or gene
        :param exp_type: str, usually is tpm or fpkm or count. default tpm
        :param group_dict: ordered dict of group info
        :param quant_method: the method to be used to quant expression.
        :param project_sn: project id
        :param task_id: task id
        :param main_id: 主表id，如果提供，则本函数不创建主表
        :param group_id: 包含分组方案信息的主表id
        :param add_distribution: 是否添加表达分布信息到数据库
        :param params: parameters dict for expression quant.
        :return: main table id
        :version:version info
        """
        # add main table info
        if main_id is None:
            name = "Gsva" + '_'
            time_now = datetime.datetime.now()
            name += time_now.strftime("%Y%m%d_%H%M%S")
            if type(params) == dict:
                params = json.dumps(params, sort_keys=True, separators=(',', ':'))
            main_info = dict(
                project_sn=project_sn,
                task_id=task_id,
                name=name,
                created_ts=datetime.datetime.now().strftime("%Y-%m-%d %H:%M:%S"),
                desc='Gsva',
                params=params,
                version="v1",
                status="start",
            )
            main_id = self.create_db_table('sg_geneset_gsva', [main_info])
        if not isinstance(main_id, ObjectId):
            main_id = ObjectId(main_id)
        df = pd.read_table(es_matrix, header=0, sep='\t')
        gsva_columns = df.columns
        columns_list = list()
        columns_list.append({'field': 'geneset', 'filter': False, 'sort':False, 'title': 'geneset', 'type': 'string'})
        for i in gsva_columns[1:]:
            columns_list.append({'field': i, 'filter': False, 'sort': False, 'title': i, 'type': 'float'})
        data_columns = {'column': columns_list, 'condition': {}}
        columns_data = json.dumps(data_columns)
        df['gsva_id'] =main_id
        df['type'] = 'heatmap'
        detail = df.to_dict('r')
        self.create_db_table('sg_geneset_gsva_detail', detail)

        main_collection = self.db['sg_geneset_gsva']
        main_collection.update({"_id": main_id},
                               {"$set": {"status": "end", 'column_data_detail': columns_data}})
        return main_id

    def add_gsva_cluster(self, cluster_output_dir, main_id=None, project_sn='medical_transcriptome',
                            task_id='medical_transcriptome',
                            params=None, title=None):
        # prepare main_table data
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
                sample_tree_list = re.findall('[(,]([^(]*?):', sample_tree)
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

        if 'gene_name' in samples:
            samples.remove('gene_name')

        if isinstance(main_id, ObjectId):
            main_id = ObjectId(main_id)
        # update main table


        tree_info_sample = dict(
            name=title,
            direction="h",
            data=sample_tree,
            type="tree",
            gsva_id=main_id,
        )
        tree_info_feature = dict(
            name=title,
            direction="v",
            data=gene_tree,
            type="tree",
            gsva_id=main_id,
        )
        tree_data = dict(name="name", condition={'type':'tree'})
        tree_data = json.dumps(tree_data)
        heatmap_data = dict(name='geneset', data=sample_tree_list, condition={'type': 'heatmap'})
        heatmap_data = json.dumps(heatmap_data)
        self.create_db_table('sg_geneset_gsva_cluster_heatmap_tree', [tree_info_sample])
        self.create_db_table('sg_geneset_gsva_cluster_heatmap_tree', [tree_info_feature])
        self.update_db_record('sg_geneset_gsva', main_id, status="end", main_id=main_id, tree_data=tree_data, heatmap_data=heatmap_data,)

    def add_gsva_diff(self, diff_output, main_id, group):
        group_dict = dict()
        with open(group, 'r') as g:
            for line in g.readlines():
                if line.startswith('#'):
                    continue
                else:
                    sample, group = line.strip().split('\t')
                    if group not in group_dict:
                        group_dict[group]=[sample]
                    else:
                        group_dict[group].append(sample)
        diff_files = glob.glob(os.path.join(diff_output, '*_vs_*.*.xls'))
        df_list = list()
        compare_list = list()
        for each in diff_files:
            diff_pd = pd.read_table(each, header=0, sep='\t')
            fname = os.path.basename(each)
            ctrl, test = re.match('(.*)_vs_(.*).limma.xls', fname).groups()
            diff_pd['compare'] = '{}|{}'.format(ctrl, test)
            compare_list.append({'field': '{}|{}'.format(ctrl, test), 'title': '{}|{}'.format(ctrl, test)})
            need_cols = ['geneset', 'log2fc', 'compare', 'fc', 'pvalue', 'padjust']
            diff_pd_new = diff_pd.loc[:,need_cols]
            df_list.append(diff_pd_new)
        df_total = pd.concat(df_list)
        df_total['gsva_id'] = main_id
        detail = df_total.to_dict('r')
        self.create_db_table('sg_geneset_gsva_diff', detail)
        gsva_columns = {'column': [{'field': 'geneset', 'filter': False, 'sort': False, 'title': 'geneset',
                                    'type': 'string'},
                                   {'field': 'log2fc', 'filter': False, 'sort': False,
                                    'title': 'log2fc',
                                    'type': 'float'},
                                   {'field': 'fc', 'filter': False, 'sort': False, 'title': 'fc',
                                    'type': 'float'},
                                   {'field': 'pvalue', 'filter': False, 'sort': False, 'title': 'pvalue',
                                    'type': 'float'},
                                   {'field': 'padjust', 'filter': False, 'sort': False, 'title': 'padjust',
                                    'type': 'float'},
                                   ], 'condition': {}}
        column_data = json.dumps(gsva_columns)
        main_collection = self.db['sg_geneset_gsva']
        main_collection.update({"_id": main_id},
                               {"$set": {"status": "end", 'column_data_diff': column_data}})
        compare_list_data = json.dumps(compare_list)
        main_collection.update({"_id": main_id},
                               {"$set": {"status": "end", 'compare': compare_list_data}})


