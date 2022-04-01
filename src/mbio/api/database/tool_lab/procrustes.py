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


class Procrustes(ApiBase):
    def __init__(self, bind_object):
        super(Procrustes, self).__init__(bind_object)
        self._project_type = 'tool_lab'

    def add_procrustes_main(self, summary, ref, query, group, main_id=None, params=None, project_sn='procrustes', task_id='procrustes'):
        # add main table info
        if main_id is None:
            name = "Procrustes" + '_'
            time_now = datetime.datetime.now()
            name += time_now.strftime("%Y%m%d_%H%M%S")
            if type(params) == dict:
                params = json.dumps(params, sort_keys=True, separators=(',', ':'))
            main_info = dict(
                project_sn=project_sn,
                task_id=task_id,
                name=name,
                created_ts=datetime.datetime.now().strftime("%Y-%m-%d %H:%M:%S"),
                desc='Procrustes',
                params=params,
                version="v1",
                status="start",
            )
            main_id = self.create_db_table('sg_procrustes', [main_info])
        if not isinstance(main_id, ObjectId):
            main_id = ObjectId(main_id)

        # detail_data required to display table
        detail_data = {
            'column': [],
            'condition': {},
        }
        page_col = ['Num included dimensions', 'Monte Carlo p-value', 'Count better', 'M2']
        data_col = ['num', 'pvalue', 'count', 'm2']
        length = len(data_col)
        for i in range(length):
            if i in [0, 2]:
                d_type = 'int'
            else:
                d_type = 'float'
            col_detail = {
                "filter": "false",
                "field": data_col[i],
                "title": page_col[i],
                "type": d_type,
                "sort": "false",
            }
            detail_data['column'].append(col_detail)

        line_data = {
            "name": "name",
            "group": 'category',
            'sub_type': 'sub_type',
            'condition': {'type': 'divide_line'}
        }

        scatter_data = {
            "name": "name",
            'data': ['x', 'y'],
            'category': 'category',
            'condition': {'type': 'scatter'}
        }

        query_dict = {
            '_id': main_id
        }
        update_dict = {
            'status': 'end',
            'main_id': main_id,
            'detail_data': json.dumps(detail_data, sort_keys=True, separators=(',', ':')),
            'divide_line_data': json.dumps(line_data, sort_keys=True, separators=(',', ':')),
            'scatter_data': json.dumps(scatter_data, sort_keys=True, separators=(',', ':')),
        }
        self.add_procrustes_graph(ref=ref, query=query, group=group, main_id=main_id)
        self.add_procrustes_detail(summary=summary, main_id=main_id)
        self.update_db_record('sg_procrustes', query_dict=query_dict, update_dict=update_dict)
        return main_id

    def add_procrustes_graph(self, ref, query, group, main_id):
        importance_list = []
        line_list = []
        site_ref = {}
        site_query = {}
        group_dict={}
        scatter_list = []
        importance = 0
        site = 0

        with open(group, 'r') as g:
            g.readline()
            for line in g:
                items = line.strip().split('\t')
                group_dict[items[0]] = items[1]

        '''
        读ref结果文件：
        获取解释度数据并写入mongo库；获取ref的样本坐标数据
        '''
        with open(ref, 'rb') as r:
            for line in r:
                arr = line.strip().split("\t")
                if not len(line) or line == "\n" or line == "\r\n":
                    importance = 0
                    site = 0
                    continue
                elif re.match("Proportion", line):
                    importance = 1
                    site = 0
                    continue
                elif importance == 1:
                    for i in range(0, len(arr)):
                        name = "PC" + str(i + 1)
                        # self.bind_object.logger.info("importance_name: %s" % name)
                        importance_data = {
                            # 'procrustes_id': ObjectId(main_id),
                            'field': str(round(float(arr[i])*100, 2))+'%',
                            'title': name,
                            # 'proportion_variance': float(arr[i]),
                            # 'type': 'importance'
                        }
                        importance_list.append(importance_data)
                    importance = 0
                elif re.match("Site", line):
                    site = 1
                    continue
                elif site == 1:
                    site_tmp = []
                    for i in range(1, len(arr)):
                        site_tmp.append(arr[i])
                    site_ref[arr[0]] = site_tmp

        query_dict = {
            '_id': main_id
        }
        update_dict = {
            'groups': json.dumps(importance_list, sort_keys=True, separators=(',', ':'))
        }
        self.update_db_record('sg_procrustes', query_dict=query_dict, update_dict=update_dict)

        '''
        读query结果文件：
        读取query的样本坐标数据
        '''
        with open(query, 'rb') as r:
            for line in r:
                arr = line.strip().split("\t")
                if not len(line) or line == "\n" or line == "\r\n":
                    site = 0
                    continue
                elif re.match("Site", line):
                    site = 1
                    continue
                elif site == 1:
                    site_tmp1 = []
                    for i in range(1, len(arr)):
                        site_tmp1.append(arr[i])
                    site_query[arr[0]] = site_tmp1
        '''
        遍历ref和query的样本坐标数据，写入mongo库
        '''
        for sample in site_ref.keys():
            site_data = {
                'procrustes_id': ObjectId(main_id),
                'name': sample,
                'type': 'divide_line',
                'category': group_dict.get(sample),
                'sub_type': 'link',
            }
            point_data_ref = {
                'procrustes_id': ObjectId(main_id),
                'name': sample,
                'type': 'scatter',
                'category': group_dict.get(sample),
            }
            point_data_query = {
                'procrustes_id': ObjectId(main_id),
                'name': sample,
                'type': 'scatter',
                'category': group_dict.get(sample),
            }
            if len(site_ref[sample]) != len(site_query[sample]):
                self.bind_object.logger.info("错误：ref和query的site坐标维度必须相同")
                self.bind_object.set_error("ERROR: site length of ref and query must be equal!")

            site_data['x1'] = site_ref[sample][0]
            site_data['x2'] = site_query[sample][0]
            site_data['y1'] = site_ref[sample][1]
            site_data['y2'] = site_query[sample][1]
            line_list.append(site_data)

            point_data_ref['x'] = site_ref[sample][0]
            point_data_ref['y'] = site_ref[sample][1]
            scatter_list.append(point_data_ref)

            point_data_query['x'] = site_query[sample][0]
            point_data_query['y'] = site_query[sample][1]
            scatter_list.append(point_data_query)

        # try:
        #     self.col_insert_data('sg_procrustes_graph', importance_list)
        # except Exception as e:
        #     self.bind_object.logger.info("导入sg_procrustes_graph:%s%s信息出错:%s" % (main_id, 'importance_list',e,))
        # else:
        #     self.bind_object.logger.info("导入sg_procrustes_graph:%s %s成功" % (main_id, 'importance_list'))

        try:
            self.col_insert_data('sg_procrustes_graph', scatter_list)
        except Exception as e:
            self.bind_object.logger.info("导入sg_procrustes_graph:%s%s信息出错:%s" % (main_id, 'scatter_list',e,))
        else:
            self.bind_object.logger.info("导入sg_procrustes_graph:%s %s成功" % (main_id, 'scatter_list'))

        try:
            self.col_insert_data('sg_procrustes_graph', line_list)
        except Exception as e:
            self.bind_object.logger.info("导入sg_procrustes_graph:%s%s信息出错:%s" % (main_id, 'line_list', e,))
        else:
            self.bind_object.logger.info("导入sg_procrustes_graph:%s %s成功" % (main_id, 'line_list'))


    def add_procrustes_detail(self, summary, main_id):
        data_list = []
        with open(summary, 'rb') as r:
            for line in r:
                if re.match("^#", line):
                    continue
                line = line.strip().split("\t")
                insert_data = {
                    'procrustes_id': ObjectId(main_id),
                    'num': int(line[2]),
                    'count': int(line[4]),
                    'm2': float(line[5]),
                    'pvalue': float(line[3])
                }
                data_list.append(insert_data)
        try:
            self.col_insert_data('sg_procrustes_detail', data_list)
        except Exception as e:
            self.bind_object.logger.info("导入sg_procrustes_detail:%s信息出错:%s" % (main_id, e,))
        else:
            self.bind_object.logger.info("导入sg_procrustes_detail:%s成功" % (main_id,))