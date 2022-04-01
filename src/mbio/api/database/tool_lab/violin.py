# !usr/bin/python
# -*- coding: utf-8 -*-
# __author__ = 'wuqin'

from collections import OrderedDict
from bson import SON
from biocluster.config import Config
import types
import re
import datetime
import pandas as pd
import os
import json
from api_base import ApiBase


class Violin(ApiBase):
    def __init__(self, bind_object):
        super(Violin, self).__init__(bind_object)
        self._db_name = 'sanger_tool_lab'

    def add_violin_detail(self, main_id, target_path):
        n_file = os.path.join(target_path + '/', 'violin.xls')
        n_table = pd.read_table(n_file, sep='\t', index_col=0)
        col = list(n_table.columns)
        for n in col:
            if n.find('group') != -1:
                col.remove(n)
        new_table = n_table[col].T.describe().T
        new_col = list(new_table.columns)
        re_table = pd.DataFrame()
        for m in new_col:
            if m.find('min') != -1:
                re_table['min'] = new_table[m]
            elif m.find('25%') != -1:
                re_table['Q1'] = new_table[m]
            elif m.find('50%') != -1:
                re_table['Q2'] = new_table[m]
            elif m.find('75%') != -1:
                re_table['Q3'] = new_table[m]
            elif m.find('max') != -1:
                re_table['max'] = new_table[m]
        re_table['irq'] = list(map(lambda x, y: x-y, re_table['Q3'], re_table['Q1']))
        # re_table.rename(columns={re_table.columns[0]: "name"}, inplace=True)
        re_table.to_csv(os.path.join(target_path, 'violin_table.txt'), sep='\t')
        new_target_file = os.path.join(target_path, 'violin_table.txt')
        with open(new_target_file, 'r') as r:
            table_datas = []
            data = r.readlines()
            table_data = data[0].strip().split('\t')
            for line in data[1:]:
                line_list = line.strip().split('\t')
                insert_data = {
                    "violin_id": self.check_objectid(main_id),
                    "type": "table",
                }
                for j in range(len(table_data)):
                    try:
                        float(line_list[j])
                    except:
                        s_data = line_list[j]
                    else:
                        s_data = float(line_list[j])
                    insert_data[table_data[j]] = s_data
                table_datas.append(insert_data)
        # 导入表格数据
        self.col_insert_data("violin_detail", table_datas)

        s_file = os.path.join(target_path + '/', 'violin_sum.xls')
        if os.path.exists(s_file):
            # 导入组内合并的小提琴画图数据
            with open(os.path.join(target_path + '/', 'violin_sum.xls'), 'r') as t1:
                s_table_datas = []
                g_lines = t1.readlines()
                g_data = g_lines[0].strip().split('\t')
                for g_line in g_lines[1:]:
                    g_line_list = g_line.strip().split('\t')
                    s_table_data = {
                        "violin_id": self.check_objectid(main_id),
                        "type": "violin",
                        "method_type": "sum"
                    }
                    s_list = []
                    for j in range(len(g_data)):
                        if j == 0:
                            s_table_data['name'] = g_line_list[j]
                            s_table_data['group'] = g_line_list[j]
                            s_table_data['x'] = g_line_list[j]
                        else:
                            try:
                                float(g_line_list[j])
                            except:
                                s_data = g_line_list[j]
                            else:
                                s_data = float(g_line_list[j])
                            s_list.append(s_data)
                        s_table_data['data'] = s_list
                    s_table_datas.append(s_table_data)
            self.col_insert_data("violin_detail", s_table_datas)
            with open(os.path.join(target_path + '/', 'violin_mean.xls'), 'r') as t2:
                m_table_datas = []
                g_lines = t2.readlines()
                g_data = g_lines[0].strip().split('\t')
                for g_line in g_lines[1:]:
                    g_line_list = g_line.strip().split('\t')
                    m_table_data = {
                        "violin_id": self.check_objectid(main_id),
                        "type": "violin",
                        "method_type": "mean"
                    }
                    m_list = []
                    for j in range(len(g_data)):
                        if j == 0:
                            m_table_data['name'] = g_line_list[j]
                            m_table_data['group'] = g_line_list[j]
                            m_table_data['x'] = g_line_list[j]
                        else:
                            try:
                                float(g_line_list[j])
                            except:
                                m_data = g_line_list[j]
                            else:
                                m_data = float(g_line_list[j])
                            m_list.append(m_data)
                        m_table_data['data'] = m_list
                    m_table_datas.append(m_table_data)
            self.col_insert_data("violin_detail", m_table_datas)
            with open(os.path.join(target_path + '/', 'violin_median.xls'), 'r') as t3:
                z_table_datas = []
                g_lines = t3.readlines()
                g_data = g_lines[0].strip().split('\t')
                for g_line in g_lines[1:]:
                    g_line_list = g_line.strip().split('\t')
                    z_table_data = {
                        "violin_id": self.check_objectid(main_id),
                        "type": "violin",
                        "method_type": "median"
                    }
                    z_list = []
                    for j in range(len(g_data)):
                        if j == 0:
                            z_table_data['name'] = g_line_list[j]
                            z_table_data['group'] = g_line_list[j]
                            z_table_data['x'] = g_line_list[j]
                        else:
                            try:
                                float(g_line_list[j])
                            except:
                                z_data = g_line_list[j]
                            else:
                                z_data = float(g_line_list[j])
                            z_list.append(z_data)
                        z_table_data['data'] = z_list
                    z_table_datas.append(z_table_data)
            self.col_insert_data("violin_detail", z_table_datas)

        # 导入组内合并为none时，或没有分组文件时小提琴画图数据
        with open(os.path.join(target_path + '/', 'violin.xls'), 'r') as v:
            v_table_datas = []
            v_lines = v.readlines()
            v_data = v_lines[0].strip().split('\t')
            for v_line in v_lines[1:]:
                v_line_list = v_line.strip().split('\t')
                v_table_data = {
                    "violin_id": self.check_objectid(main_id),
                    "type": "violin",
                    "name": v_line_list[0],
                    "group": v_line_list[len(v_data) - 1],
                    "x": v_line_list[0]
                }
                v_list = []
                if v_data[-1] == "group":
                    for j in range(1, len(v_data)-1):
                        try:
                            float(v_line_list[j])
                        except:
                            f_data = v_line_list[j]
                        else:
                            f_data = float(v_line_list[j])
                        v_list.append(f_data)
                        v_table_data['data'] = v_list
                        if os.path.exists(s_file):
                            v_table_data['method_type'] = "none"
                else:
                    for j in range(1, len(v_data)):
                        try:
                            float(v_line_list[j])
                        except:
                            f_data = v_line_list[j]
                        else:
                            f_data = float(v_line_list[j])
                        v_list.append(f_data)
                        v_table_data['data'] = v_list
                v_table_datas.append(v_table_data)
        self.col_insert_data("violin_detail", v_table_datas)

        # 更新主表
        column = [
            {"field": "name", "type": "string", "sort": "false", "title": "id"},
            {"field": "min", "type": "float", "sort": "false", "title": "min"},
            {"field": "Q1", "type": "float", "sort": "false", "title": "Q1"},
            {"field": "Q2", "type": "float", "sort": "false", "title": "Q2"},
            {"field": "Q3", "type": "float", "sort": "false", "title": "Q3"},
            {"field": "max", "type": "float", "sort": "false", "title": "max"},
            {"field": "irq", "type": "float", "sort": "false", "title": "IRQ"}
        ]
        table_data_update = {'column': column, 'condition': {'type': 'table'}}
        update_info = {
            "main_id": self.check_objectid(main_id),
            "table_data": json.dumps(table_data_update),
        }
        if os.path.exists(s_file):
            violin_data_update = {"name": "name", "data": "data", "condition": {"type": "violin"},
                                  "data_option": ["none", "sum", "mean", "median"]}
        else:
            violin_data_update = {"name": "name", "data": "data", "condition": {"type": "violin"}}
        update_info["violin_data"] = json.dumps(violin_data_update)
        self.update_db_record("violin", {"_id": self.check_objectid(main_id)}, update_info)


if __name__ == "__main__":
    a = Violin(None)
    a.add_violin_detail("5ea7b87917b2bf0f61841111", "/mnt/ilustre/users/sanger-dev/workspace/20200610/Violin_tsg_3421_0610161957838831_1046/output/violin")
    # a.add_violin_detail("5ea7b87917b2bf0f61841111", "/mnt/ilustre/users/sanger-dev/workspace/20200610/Violin_tsg_3421_0610162903961746_6528/output/violin")
