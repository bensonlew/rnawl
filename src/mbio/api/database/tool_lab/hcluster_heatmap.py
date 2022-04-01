# -*- coding: utf-8 -*-
# __author__ = 'zoujiaxun'
from bson.objectid import ObjectId
import datetime
import json
import pickle
import unittest
import pandas as pd
from mbio.api.database.whole_transcriptome.api_base import ApiBase
import re
import os


class HclusterHeatmap(ApiBase):
    def __init__(self, bind_object):
        super(HclusterHeatmap, self).__init__(bind_object)
        self._project_type = 'tool_lab'

    def add_hcluster_heatmap(self, sampletree, featuretree, heatmap, group, project_sn='tool_lab', main_id=None,
                     task_id='tool_lab', title=None):
        group_dict = dict()
        group_list = list()
        gene_cluster, sample_cluster = False, False
        try:
            with open(group, 'r') as groups:
                for line in groups:
                    if line.startswith("#"):
                        continue
                    sample, group = line.strip().split('\t')
                    if group not in group_dict:
                        group_dict[group] = [sample]
                    else:
                        group_dict[group].append(sample)
            for group in group_dict:
                group_list.append({'groupname': group, 'value':group_dict[group]})
        except:
            pass
        if main_id is None:
            name = "Hcluster_heatmap"+'_'
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
                desc='Cluster heatmap tree main table',
                params=[],
                status="start",
            )
            main_id = self.create_db_table('cluster_heatmap', [main_info])
        else:
            main_id = ObjectId(main_id)

        if type(main_id) == str or type(main_id) == bytes or type(main_id) == unicode:
            main_id = ObjectId(main_id)
        if os.path.isfile(sampletree):
            target_file_sample = sampletree
            sample_cluster = True
        if os.path.isfile(featuretree):
            target_file_feature = featuretree
            gene_cluster = True

        if sample_cluster:
            with open(target_file_sample) as f:
                sample_tree = f.readline().strip()
                sample_tree_list = re.findall('[(,]([^(]*?):', sample_tree)
                samples = f.readline().strip().split(";")

        if gene_cluster:
            with open(target_file_feature) as f:
                feature_tree = f.readline().strip()
                feature_tree_list = re.findall('[(,]([^(]*?):', feature_tree)
                samples = f.readline().strip().split(";")

        heatmap = heatmap
        df = pd.read_table(heatmap, index_col=0, sep='\t')
        sample_list = df.columns.tolist()
        if gene_cluster:
            df1 = df.reindex(feature_tree_list)
        else:
            df1 = df
        df2 = df1.reset_index()
        df2['cluster_id'] = main_id
        df2['type'] = 'heatmap'
        self.create_db_table('cluster_heatmap_detail', df2.to_dict('r'))

        # add detail info
        if gene_cluster and sample_cluster:
            if group_list is None:
                tree_info_sample = dict(
                    name=title,
                    direction="h",
                    data=sample_tree,
                    # data2=feature_tree,
                    type="tree",
                    cluster_id=main_id,
                )
            else:
                tree_info_sample = dict(
                    name=title,
                    direction="h",
                    data=sample_tree,
                    # data2=feature_tree,
                    type="tree",
                    cluster_id=main_id,
                    group=group_list
                )
            tree_info_feature = dict(
                name=title,
                direction="v",
                data=feature_tree,
                # data2=feature_tree,
                type="tree",
                cluster_id=main_id,
                # group_dict=group_dict
            )
        tree_data = dict(name="name", condition={'type':'tree'})
        tree_data = json.dumps(tree_data)

        # sample_list = list()
        # for group in group_dict:
        #     sample_list.extend(group_dict[group])
        # heatmap_data = dict(name='feature', data=sample_list, condition={'type':'heatmap'})
        if sample_cluster:
            heatmap_data = dict(name='feature', data=sample_tree_list, condition={'type':'heatmap'})
        else:
            heatmap_data = dict(name='feature', data=sample_list, condition={'type': 'heatmap'})
        heatmap_data = json.dumps(heatmap_data)
        # sample_bar = re.findall('[(,]([^(]*?):', sample_tree)
        # print sample_bar
        if gene_cluster and sample_cluster:
            heatmap_bar_list = list()
            # for s in sample_bar:
            #     heatmap_bar_list.extend(filter(lambda k: s in group_dict[k], group_dict))
            # print heatmap_bar_list

            for k in group_dict:
                for j in group_dict[k]:
                    heatmap_bar = dict(
                        name = j,
                        direction = 'h',
                        type = 'heatmap_bar',
                        level_color = k,
                        cluster_id = main_id
                    )
                    heatmap_bar_list.append(heatmap_bar)
        heatmap_bar_data = dict(name = 'name', condition = {'type': 'heatmap_bar'})
        heatmap_bar_data = json.dumps(heatmap_bar_data)
        if gene_cluster and sample_cluster:
            self.create_db_table('cluster_heatmap_tree', [tree_info_sample])
            self.create_db_table('cluster_heatmap_tree', [tree_info_feature])
            self.create_db_table('cluster_heatmap_bar', heatmap_bar_list)
        self.update_db_record('cluster_heatmap', main_id, status="end", main_id=main_id, tree_data=tree_data, heatmap_data=heatmap_data, heatmap_bar_data=heatmap_bar_data)
        # else:
        #     self.update_db_record('cluster_heatmap', main_id, status="end", main_id=main_id, heatmap_data=heatmap_data)
        return main_id
