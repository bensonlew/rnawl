# -*- coding: utf-8 -*-
# __author__ = 'zoujiaxun'
from bson.objectid import ObjectId
import datetime
import json
import pickle
import unittest
import pandas as pd
from mbio.api.database.whole_transcriptome.api_base import ApiBase


class Venn(ApiBase):
    def __init__(self, bind_object):
        super(Venn, self).__init__(bind_object)
        self._project_type = 'tool_lab'

    def add_venn(self, venn_graph, venn_table=None, project_sn='tool_lab', exp_level='T', main_id=None,
                     quant_method='RSEM', task_id='tool_lab', params=None):
        if main_id is None:
            name = "Venn"+'_'
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
                desc='venn main table',
                params=[],
                status="start",
            )
            main_id = self.create_db_table('venn', [main_info])
        else:
            main_id = ObjectId(main_id)

        if type(main_id) == str or type(main_id) == bytes or type(main_id) == unicode:
            main_id = ObjectId(main_id)
        graph_pd = pd.read_table(venn_graph, header=0, sep='\t', keep_default_na=False)
        graph_pd.columns = ["name", "seqs"]
        detail_dict_list = graph_pd.to_dict('records')
        # category = graph_pd['name'].tolist()
        # if venn_table:
        #     table_pd = pd.read_table(venn_table, header=None, sep='\t', keep_default_na=False)
        #     table_pd.columns = ["combination", "num", "only_list"]
        #     detail_dict_list += table_pd.to_dict('records')
        self.create_db_table('venn_detail', detail_dict_list, tag_dict={'venn_id': main_id, 'type': 'venn'})
        venn_data = dict(names='seqs',category='name', condition={'type':'venn'})
        venn_data = json.dumps(venn_data,  sort_keys=True, separators=(',', ':'))
        self.update_db_record('venn', main_id, status="end", main_id=main_id, venn_data=venn_data)
        return main_id
