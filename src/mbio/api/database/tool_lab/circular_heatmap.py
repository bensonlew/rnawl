# -*- coding: utf-8 -*-
# __author__ = 'zhangyitong'

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
from api_base import ApiBase


class CircularHeatmap(ApiBase):
    def __init__(self, bind_object):
        super(CircularHeatmap, self).__init__(bind_object)
        self._project_type = 'tool_lab'

    def add_circular_heatmap(self, sampletree, featuretree, heatmap, group, project_sn='tool_lab', main_id=None,
                             task_id='tool_lab', title=None):
        group_dict = dict()
        group_list = list()
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
                group_list.append({'groupname': group, 'value': group_dict[group]})
        except:
            pass
        if main_id is None:
            name = "Circular_Heatmap" + '_'
            time_now = datetime.datetime.now()
            name += time_now.strftime("%Y%m%d_%H%M%S")
            # if type(params) == dict:
            #     params = json.dumps(params, sort_keys=True, separators=(',', ':'))
            main_info = dict(
                project_sn=project_sn,
                task_id=task_id,
                version="v1",
                name=name,
                created_ts=time_now.strftime('%Y-%m-%d %H:%M:%S'),
                desc='circular heatmap main table',
                params=[],
                status="start",
            )
            main_id = self.col_insert_data('sg_circular_heatmap', [main_info])
        else:
            main_id = ObjectId(main_id)

        if sampletree:
            with open(sampletree) as f:
                sample_tree = f.readline().strip()
                sample_tree_list = re.findall('[(,]([^(]*?):', sample_tree)
                samples = f.readline().strip().split(";")

        if featuretree:
            with open(featuretree) as f:
                feature_tree = f.readline().strip()
                feature_tree_list = re.findall('[(,]([^(]*?):', feature_tree)
                # features = f.readline().strip().split(";")

        heatmap = heatmap
        df = pd.read_table(heatmap, index_col=0, sep='\t')
        df.drop(columns=['metab_id'], inplace=True)
        sample_list = df.columns.tolist()
        if featuretree:
            df1 = df.reindex(feature_tree_list)
        else:
            df1 = df
        df2 = df1.reset_index()
        df2['cluster_id'] = main_id
        df2['type'] = 'table'
        df2['Metabolite'] = df2.apply(lambda x: x['Metabolite'].replace('_@', '.'), axis=1)
        df2.rename(columns={'Metabolite': 'metabolite'}, inplace=True)
        self.col_insert_data('sg_circular_heatmap_table', df2.to_dict('r'))

        df_t = df.T
        metab_list = df_t.columns.tolist()
        if sampletree:
            df1 = df_t.reindex(sample_tree_list)
        else:
            df1 = df_t
        df2 = df1.reset_index()
        df2['cluster_id'] = main_id
        df2['type'] = 'heatmap'
        df2.rename(columns={'index': 'sample'}, inplace=True)
        self.col_insert_data('sg_circular_heatmap_detail', df2.to_dict('r'))

        # add detail info
        if sampletree:
            if group_list is None:
                tree_info_sample = dict(
                    name=title,
                    direction="v",
                    data=sample_tree,
                    type="tree",
                    cluster_id=main_id,
                )
            else:
                tree_info_sample = dict(
                    name=title,
                    direction="v",
                    data=sample_tree,
                    type="tree",
                    cluster_id=main_id,
                    group=group_list
                )
            self.col_insert_data('sg_circular_heatmap_tree', [tree_info_sample])
        if featuretree:
            tree_info_feature = dict(
                name=title,
                direction="circle",
                data=feature_tree,
                type="tree",
                cluster_id=main_id,
            )
            self.col_insert_data('sg_circular_heatmap_tree', [tree_info_feature])
        tree_data = dict(name="name", condition={'type': 'tree'})
        tree_data = json.dumps(tree_data)

        if featuretree:
            heatmap_data = dict(name='sample', data=feature_tree_list, condition={'type': 'heatmap'})
        else:
            heatmap_data = dict(name='sample', data=metab_list, condition={'type': 'heatmap'})
        heatmap_data = json.dumps(heatmap_data)

        # if featuretree and sampletree:
        #     heatmap_bar_list = list()
        #
        #     for k in group_dict:
        #         for j in group_dict[k]:
        #             heatmap_bar = dict(
        #                 name=j,
        #                 direction='circle',
        #                 type='heatmap_bar',
        #                 level_color=k,
        #                 cluster_id=main_id
        #             )
        #             heatmap_bar_list.append(heatmap_bar)
        # heatmap_bar_data = dict(name='name', condition={'type': 'heatmap_bar'})
        # heatmap_bar_data = json.dumps(heatmap_bar_data)
        # self.col_insert_data('sg_circular_heatmap_bar', heatmap_bar_list)

        detail_data = {
            'column': [],
            'condition': {'type': 'table'},
        }
        data_col = sample_list
        length = len(data_col) + 1
        for i in range(length):
            if i == 0:
                d_type = 'string'
                title = 'Row name'
                field = 'metabolite'
            else:
                d_type = 'float'
                title = data_col[i-1]
                field = data_col[i-1]
            col_detail = {
                "filter": "false",
                "field": field,
                "title": title,
                "type": d_type,
                "sort": "false",
            }
            detail_data['column'].append(col_detail)

        query_dict = {
            '_id': main_id
        }
        update_dict = {
            'status': 'end',
            'main_id': main_id,
            'detail_data': json.dumps(detail_data, sort_keys=True, separators=(',', ':')),
            # 'heatmap_bar_data': heatmap_bar_data,
            'heatmap_data': heatmap_data,
            'tree_data': tree_data,
        }

        self.update_db_record('sg_circular_heatmap', query_dict=query_dict, update_dict=update_dict)
        return main_id