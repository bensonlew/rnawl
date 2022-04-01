# -*- coding: utf-8 -*-
# __author__ = 'zzg'

from bson.objectid import ObjectId
import datetime
import json
import os
import re
import pickle
import unittest
import pandas as pd
from mbio.api.database.whole_transcriptome.api_base import ApiBase


class CorHeatmap(ApiBase):
    def __init__(self, bind_object):
        super(CorHeatmap, self).__init__(bind_object)
        self._project_type = 'tool_lab'

    def add_corheatmap(self, CorHeatmap_path, project_sn='tool_lab', main_id=None, task_id='tool_lab', params=None):
        if main_id is None:
            name = "corheatmap"+'_'
            time_now = datetime.datetime.now()
            name += time_now.strftime("%Y%m%d_%H%M%S")
            main_info = dict(
                project_sn=project_sn,
                task_id=task_id,
                version="v2",
                name=name,
                created_ts=time_now.strftime('%Y-%m-%d %H:%M:%S'),
                desc='corheatmap',
                params= params if params else "null",
                status="start",
            )
            main_id = self.create_db_table('sg_cor_heatmap', [main_info])
        else:
                main_id = ObjectId(main_id)

        if type(main_id) == str or type(main_id) == bytes or type(main_id) == unicode:
            main_id = ObjectId(main_id)

        if os.path.exists(CorHeatmap_path + "/correlation_spearmanr.xls"):
            correction_file = CorHeatmap_path + "/correlation_spearmanr.xls"
            pvalue_file = CorHeatmap_path + "/pvalue_spearmanr.xls"
        else:
            correction_file = CorHeatmap_path + "/correlation_pearsonr.xls"
            pvalue_file = CorHeatmap_path + "/pvalue_pearsonr.xls"

        if os.path.exists(CorHeatmap_path + "/row_tree.tre"):
            with open(CorHeatmap_path + "/row_tree.tre") as f:
                row_tree = f.readline().strip()
                row_tree_list = re.findall('[(,]([^(]*?):', row_tree)
                tree_info_row = dict(
                    name="row_tree",
                    direction="h",
                    data=row_tree,
                    type="tree",
                    cor_id=main_id,
                )
                self.create_db_table('sg_cor_heatmap_detail', [tree_info_row])
        if os.path.exists(CorHeatmap_path + "/column_tree.tre"):
            with open(CorHeatmap_path + "/column_tree.tre") as f:
                column_tree = f.readline().strip()
                column_tree_list = re.findall('[(,]([^(]*?):', column_tree)
                tree_info_column = dict(
                    name="column_tree",
                    direction="v",
                    data=column_tree,
                    type="tree",
                    cor_id=main_id,
                )
                self.create_db_table('sg_cor_heatmap_detail', [tree_info_column])
        """
         column_tree_new_name = []
        column_tree_name_dict = {}
        for i in range(len(column_tree_list)):
            column_tree_new_name.append("name" + str(i + 1))
            column_tree_name_dict["name" + str(i + 1)] = column_tree_list[i]
        column_name = json.dumps(column_tree_name_dict)
        """
        column_name = []
        with open(correction_file,"r") as h:
            a = h.readlines()
            for i in range(len(a[1:])):
                column_name.append(a[i+1].strip("\n").split("\t")[0])
            df1 = pd.read_table(correction_file, index_col=0, sep='\t')
            df_cor1 = df1.reindex(column_name)
            #df_cor.index =  column_name
            df_cor = df_cor1.reset_index()
            df_cor['cor_id'] = main_id
            df_cor['type'] = 'heatmap'
            #df_cor['name'] = column_tree_new_name
            self.create_db_table('sg_cor_heatmap_detail', df_cor.to_dict('r'))
            df2 = pd.read_table(pvalue_file, index_col=0, sep='\t')
            df_pva1 = df2.reindex(column_name)
            #df_pva.index = column_tree_new_name
            df_pva = df_pva1.reset_index()
            df_pva['cor_id'] = main_id
            df_pva['type'] = 'heatmap_asterisk'
            #df_pva['name'] = column_tree_new_name
            self.create_db_table('sg_cor_heatmap_detail', df_pva.to_dict('r'))
        '''
        heatmap_data = {"heatmap_data": [{"name": "name", "data": row_name}],
                        "condition": {"type": "heatmap"}}
        heatmap_asterisk_data = {"heatmap_asterisk_data": [{"name": "sample_name", "data": row_name}],
                                 "condition": {"type": "heatmap_asterisk"}}
        tree_data = {"tree_data": [{"name": "tree_name"}],
                     "condition": {"type": "tree "}}
        heatmap_data_info = json.dumps(heatmap_data, sort_keys=False, separators=(',', ':'))
        heatmap_asterisk_data_info = json.dumps(heatmap_asterisk_data, sort_keys=False, separators=(',', ':'))
        tree_data_info = json.dumps(tree_data, sort_keys=False, separators=(',', ':'))
        '''
        with open(correction_file,"r") as h:
            a = h.readlines()
            row_name = a[0].strip("\n").split("\t")[1:]
            row_name.insert(0,"name")
            heatmap_data = dict(name='name', data=row_name, condition={'type': 'heatmap'})
            heatmap_data = json.dumps(heatmap_data)
            heatmap_asterisk_data =  dict(name='name', data=row_name, condition={'type': 'heatmap_asterisk'})
            heatmap_asterisk_data = json.dumps(heatmap_asterisk_data)
            tree_data = dict(name='tree_name',condition={'type': 'tree'})
            tree_data = json.dumps(tree_data)
            self.update_db_record('sg_cor_heatmap', main_id, status="end", main_id=main_id, heatmap_data=heatmap_data, heatmap_asterisk_data=heatmap_asterisk_data, tree_data=tree_data,) #column_name=column_tree_name_dict,
        return main_id
