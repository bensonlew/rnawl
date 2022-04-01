# !usr/bin/python
# -*- coding: utf-8 -*-
from collections import OrderedDict
from bson import SON
from biocluster.config import Config
import types
import re
import datetime
import os
import json
from api_base import ApiBase


class Correlation(ApiBase):
    """..."""
    def __init__(self, bind_object):
        super(Correlation, self).__init__(bind_object)

    def insert_detail_table(self, main_id, heatmap_file, tree_file, tree_2_file, pvalue_file):
        """
        这里着重说明下，表格是如何导入的
        table_data： {"column":
         [{"field": "name", "title": "name", "filter": "false", "sort": "false", "type": "string"},
         {"field": "y1", "title": "y1", "filter": "false", "sort": "false", "type": "float"}],
         "condition": {'type': "heatamap"}}
         column: 字段信息
         filed：字段名
         title:页面展示的名称
         filter：是否要过滤 （true or false）
         sort：是否能排序（true or false）
         type：数值类型（string， float， int）
        :param main_id:
        :param heatmap_file:
        :param tree_file: corr_col 垂直树
        :param tree_2_file: corr_row  水平树
        :return:
        """
        insert_datas = []
        column = [{"field": "name", "title": "name", "filter": "false", "sort": "false", "type": "string"}]
        p_column = [{"field": "name", "title": "name", "filter": "false", "sort": "false", "type": "string"}]
        with open(heatmap_file, 'r') as r:
            data = r.readlines()
            heatmap_data = data[0].strip().split('\t')
            for j in range(len(heatmap_data)):
                oneinfo = {"field": heatmap_data[j], "title": heatmap_data[j], "filter": "false",
                           "sort": "false", "type": "float"}
                column.append(oneinfo)
                insert_data = {
                    'name': heatmap_data[j],
                    "correlation_id": self.check_objectid(main_id),
                    "type": "heatmap"
                }
                for line in data[1:]:
                    temp = line.strip().split('\t')
                    insert_data[temp[0]] = float(temp[j + 1])
                insert_datas.append(insert_data)
        with open(pvalue_file, 'r') as p:
            p_data = p.readlines()
            pvalue_data = p_data[0].strip().split('\t')
            for i in range(len(pvalue_data)):
                p_oneinfo = {"field": pvalue_data[i], "title": pvalue_data[i], "filter": "false",
                           "sort": "false", "type": "float"}
                p_column.append(p_oneinfo)
                p_insert_data = {
                    'name': pvalue_data[i],
                    "correlation_id": self.check_objectid(main_id),
                    "type": "heatmap_asterisk"
                }
                for p_line in p_data[1:]:
                    p_temp = p_line.strip().split('\t')
                    p_insert_data[p_temp[0]] = float(p_temp[i + 1])
                insert_datas.append(p_insert_data)
        if tree_file:
            with open(tree_file, 'r') as r1:
                data = r1.readlines()[0].strip()
                insert_datas.append({
                    "name": 'tree_col',
                    "type": "tree",
                    "direction": "v",
                    "correlation_id": self.check_objectid(main_id),
                    "data": data
                })
        if tree_2_file:
            with open(tree_2_file, 'r') as r2:
                data = r2.readlines()[0].strip()
                insert_datas.append({
                    "name": 'tree_row',
                    "type": "tree",
                    "direction": "h",
                    "correlation_id": self.check_objectid(main_id),
                    "data": data
                })
        self.col_insert_data("sg_correlation_detail", insert_datas)
        heatmap_data_insert = {
            "name": "name",
            "data": heatmap_data,
            "condition": {'type': "heatmap"}
        }
        tree_data_insert = {
            "name": "name",
            "data": ['data'],
            "condition": {'type': "tree"}
        }
        table_data_insert = {
            "column": column,
            "condition": {'type': "heatmap"}
        }
        pvalue_data_insert = {
            "column": p_column,
            "condition": {'type': "heatmap_asterisk"}
        }
        asterisk_data_insert = {
            "name": "name",
            "data": pvalue_data,
            "condition": {'type': "heatmap_asterisk"}
        }
        update_dict = {
            "heatmap_data": json.dumps(heatmap_data_insert),
            "tree_data": json.dumps(tree_data_insert),
            "table_data": json.dumps(table_data_insert),
            "pvalue_data": json.dumps(pvalue_data_insert),
            "asterisk_data": json.dumps(asterisk_data_insert)
        }
        self.update_db_record("sg_correlation", {"_id": self.check_objectid(main_id)}, update_dict)
        pass


if __name__ == "__main__":
    a = Correlation(None)
    a.insert_detail_table('5ea91cbf17b2bf65274c6666', '/mnt/ilustre/users/sanger-dev/workspace/20210322/Correlation_untl_cqslolqk37e11egl68akl2_0322090421886796_3120/output/correlation/correlation_matrix.xls', 
        '/mnt/ilustre/users/sanger-dev/workspace/20210322/Correlation_untl_cqslolqk37e11egl68akl2_0322090421886796_3120/output/correlation/corr_col.tre',
        '/mnt/ilustre/users/sanger-dev/workspace/20210322/Correlation_untl_cqslolqk37e11egl68akl2_0322090421886796_3120/output/correlation/corr_row.tre',
        '/mnt/ilustre/users/sanger-dev/workspace/20210322/Correlation_untl_cqslolqk37e11egl68akl2_0322090421886796_3120/output/correlation/pvalue_matrix.xls')
    # a.table_in("/mnt/ilustre/users/sanger-dev/workspace/20190827/Single_correlation20190827092037/Correlation/output/correlation_matrix.xls")
    # a.insert_table("/mnt/ilustre/users/sanger-dev/workspace/20190827/Single_correlation20190827092037/Correlation/output/correlation_matrix.xls",
    #                corr_result, corr)
    # a.correlation_in("/mnt/ilustre/users/sanger-dev/workspace/20190827/Single_correlation20190827092037/Correlation/output/correlation_matrix.xls",
    #                  "/mnt/ilustre/users/sanger-dev/workspace/20190827/Single_correlation20190827092037/Correlation/output/cluster.tre")