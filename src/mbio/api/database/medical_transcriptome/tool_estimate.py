# -*- coding: utf-8 -*-
# __author__ = 'zoujiaxun'
# last_modify:20201012
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
from mbio.api.database.medical_transcriptome.api_base import ApiBase

class ToolEstimate(ApiBase):
    def __init__(self, bind_object):
        super(ToolEstimate, self).__init__(bind_object)

    def add_estimate(self, main_id, estimate_score, params=None, project_sn='medical_transcriptome', task_id='medical_transcriptome'):
        # add main table info
        if main_id is None:
            name = "Estimate" + '_'
            time_now = datetime.datetime.now()
            name += time_now.strftime("%Y%m%d_%H%M%S")
            if type(params) == dict:
                params = json.dumps(params, sort_keys=True, separators=(',', ':'))
            main_info = dict(
                project_sn=project_sn,
                task_id=task_id,
                name=name,
                created_ts=datetime.datetime.now().strftime("%Y-%m-%d %H:%M:%S"),
                desc='Estimate',
                params=params,
                version="v1",
                status="start",
            )
            main_id = self.create_db_table('sg_tool_estimate', [main_info])
        if not isinstance(main_id, ObjectId):
            main_id = ObjectId(main_id)
        df = pd.read_table(estimate_score, header=0, sep='\t', skiprows=2)
        df['estimate_id'] = main_id
        detail = df.to_dict('r')
        self.create_db_table('sg_tool_estimate_detail', detail)
        self.update_db_record('sg_tool_estimate', main_id, status='end', main_id=main_id)
        return main_id

    def add_estimate_cluster(self, cluster_output_dir, main_id=None, project_sn='medical_transcriptome',
                            task_id='medical_transcriptome',
                            params=None):
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
        self.update_db_record('sg_tool_estimate', main_id,
                              sample_cluster=sample_cluster)

        tree_info = dict(
            sample_tree=sample_tree,
            estimate_id=main_id,
            samples=samples,
        )

        self.create_db_table('sg_tool_estimate_cluster_tree', [tree_info])
