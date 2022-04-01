# -*- coding: utf-8 -*-
# __author__ = 'zoujiaxun'
from bson.objectid import ObjectId
import datetime
import json
import pickle
import unittest
import pandas as pd
from mbio.api.database.whole_transcriptome.api_base import ApiBase


class Hcluster(ApiBase):
    def __init__(self, bind_object):
        super(Hcluster, self).__init__(bind_object)
        self._project_type = 'tool_lab'

    def add_hcluster(self, tree, group, project_sn='tool_lab', exp_level='T', main_id=None,
                     quant_method='RSEM', task_id='tool_lab', params=None, main_title=None):
        group_dict = dict()
        with open(group, 'r') as groups:
            for line in groups:
                if line.startswith("#"):
                    continue
                sample, group = line.strip().split('\t')
                if group not in group_dict:
                    group_dict[group] = [sample]
                else:
                    group_dict[group].append(sample)
        group_list = list()
        for group in group_dict:
            group_list.append({'groupname': group, 'value':group_dict[group]})

        if main_id is None:
            name = "Hcluster"+'_'
            time_now = datetime.datetime.now()
            name += time_now.strftime("%Y%m%d_%H%M%S")
            # if type(params) == dict:
            #     params = json.dumps(params, sort_keys=True, separators=(',', ':'))
            main_info = dict(
                project_sn=project_sn,
                task_id=task_id,
                version="v3",
                name=name,
                created_ts=time_now.strftime('%Y-%m-%d %H:%M:%S'),
                desc='Cluster tree main table',
                params=[],
                status="start",
            )
            main_id = self.create_db_table('cluster', [main_info])
        else:
            main_id = ObjectId(main_id)

        if type(main_id) == str or type(main_id) == bytes or type(main_id) == unicode:
            main_id = ObjectId(main_id)
        target_file = tree
        with open(target_file) as f:
            sample_cluster = True
            sample_tree = f.readline().strip()
            samples = f.readline().strip().split(";")
        # add detail info
        tree_info = dict(
            name=main_title,
            direction="h",
            data=sample_tree,
            type="tree",
            cluster_id=main_id,
            group=group_list
        )
        tree_data = dict(name="name", condition={'type':'tree'})
        tree_data = json.dumps(tree_data)
        self.create_db_table('cluster_detail', [tree_info])
        self.update_db_record('cluster', main_id, status="end", main_id=main_id, tree_data=tree_data)
        return main_id
