# !/usr/bin/python
# -*- coding: utf-8 -*-
from bson.objectid import ObjectId
# import types
# import re
# import os
import json
import pandas as pd
# import numpy as np
# from scipy import stats
# import math
# from sklearn import decomposition
# import scipy.cluster.hierarchy as sch
# import fastcluster as hclust
# from collections import OrderedDict
# import unittest
import datetime
# import glob
from api_base import ApiBase


class AllExp(ApiBase):
    def __init__(self, bind_object):
        super(AllExp, self).__init__(bind_object)

    def add_exp_venn(self, venn_graph, venn_table=None, project_sn='ref_rna', exp_level='T', main_id=None,
                     quant_method='RSEM', task_id='ref_rna', params=None):
        """
        add venn analysis info
        :param venn_graph: venn_graph.xls resulted from express_venn tool
        :param venn_table: venn_table.xls resulted from express_venn tool
        :param quant_method: exp quant method
        :param project_sn: project id
        :param task_id: task id
        :param main_id: 主表id，如果提供，则本函数不创建主表
        :param exp_level: T or G
        :param params: string of parameters dict, designed for judgement of whether the task is repeated.
        :return: main table id
        """
        # add main table info
        if main_id is None:
            name = "ExpVenn" + '_' + exp_level + '_' + quant_method + '_'
            time_now = datetime.datetime.now()
            name += time_now.strftime("%Y%m%d_%H%M%S")
            if type(params) == dict:
                params = json.dumps(params, sort_keys=True, separators=(',', ':'))
            main_info = dict(project_sn=project_sn, task_id=task_id, name=name,
                created_ts=time_now.strftime('%Y-%m-%d %H:%M:%S'), desc='venn main table', params=params, status="start")
            main_id = self.create_db_table('sg_exp_venn', [main_info])
        else:
            main_id = ObjectId(main_id)
        # add detail table info
        graph_pd = pd.read_table(venn_graph, header=0, sep='\t')
        graph_pd.columns = ["name", "seqs"]
        detail_dict_list = graph_pd.to_dict('records')
        if venn_table:
            table_pd = pd.read_table(venn_table, header=None, sep='\t')
            table_pd.columns = ["combination", "num", "only_list"]
            detail_dict_list += table_pd.to_dict('records')
        self.create_db_table('sg_exp_venn_detail', detail_dict_list, tag_dict={'venn_id': main_id})
        self.update_db_record('sg_exp_venn', main_id, status="end", main_id=main_id)
        return main_id

