# -*- coding: utf-8 -*-
# __author__ = 'wuqin'

from collections import OrderedDict
from bson import SON
from biocluster.config import Config
import pandas as pd
from io import StringIO
import types
import re
import datetime
import os
import json
from api_base import ApiBase


class Pca(ApiBase):
    def __init__(self, bind_object):
        super(Pca, self).__init__(bind_object)
        self._db_name = 'sanger_tool_lab'

    def add_pca_detail(self, main_id, output_dir):
        """
        导入pca画散点图数据
        """
        target_file = os.path.join(output_dir, 'pca_sites.xls')
        if os.path.exists(os.path.join(output_dir, 'group.xls')):
            group_file = os.path.join(output_dir, 'group.xls')
            site = pd.read_table(target_file)
            group = pd.read_table(group_file)
            site.rename(columns={site.columns[0]: "name"}, inplace=True)
            group.rename(columns={group.columns[0]: "name", group.columns[1]: "group"}, inplace=True)
            g_key = group.columns[0]
            new_site = pd.merge(site, group, on=g_key)
            new_site.to_csv(os.path.join(output_dir, 'pca_sites_group.xls'), sep='\t', index=0)
            new_target_file = os.path.join(output_dir, 'pca_sites_group.xls')
        else:
            site = pd.read_table(target_file)
            site.rename(columns={site.columns[0]: "name"}, inplace=True)
            site.to_csv(os.path.join(output_dir, 'pca_sites_new.xls'), sep='\t', index=0)
            new_target_file = os.path.join(output_dir, 'pca_sites_new.xls')
        with open(new_target_file, 'r') as r:
            insert_datas = []
            data = r.readlines()
            pca_data = data[0].strip().split('\t')
            for line in data[1:]:
                line_list = line.strip().split('\t')
                insert_data = {
                    "pca_id": self.check_objectid(main_id),
                    "type": "scatter",
                }
                for j in range(len(pca_data)):
                    try:
                        float(line_list[j])
                    except:
                        s_data = line_list[j]
                    else:
                        if float(line_list[j]) < 0.000001:
                            s_data = '{:e}'.format(float(line_list[j]))
                        else:
                            s_data = float(line_list[j])
                    insert_data[pca_data[j]] = s_data
                    if 'group' not in pca_data:
                        insert_data['group'] = ''
                insert_datas.append(insert_data)
        self.col_insert_data("pca_detail", insert_datas)
        new_pca_data = pca_data
        if 'name' in new_pca_data:
            new_pca_data.remove('name')
        if 'group' in new_pca_data:
            new_pca_data.remove('group')
        pca_data_insert = {
            "name": "name",
            "data": new_pca_data,
            "category": "group",
            "condition": {'type': "scatter"}
        }
        update_dict = {
            "scatter_data": json.dumps(pca_data_insert),
        }
        self.update_db_record("pca", {"_id": self.check_objectid(main_id)}, update_dict)

    def add_pca_ellipse_detail(self, main_id, output_dir):
        ellipse = os.path.join(output_dir, 'ellipse.xls')
        with open(ellipse, 'r') as s:
            s_data = s.readlines()
            s_line_list0 = s_data[0].strip().split('\t')
            s_insert_datas1 = []
            if s_line_list0[1] == "All":
                # 没有分组时
                for s_line in s_data[1:]:
                    s_line_list = s_line.strip().split('\t')
                    s_insert_data1 = {
                        "pca_id": self.check_objectid(main_id),
                        "type": "ellipse",
                        "name": s_line_list0[1],
                        "group": s_line_list0[1],
                        "method": "sample"
                    }
                    data_list = s_line_list[1].split(",")
                    if len(data_list) == 6:
                        # 数据出来为 m1, c11, c12, m2, c21, c22
                        s_data1 = str(float(data_list[0])) + ',' + str(float(data_list[1])) + ',' + \
                                  str(float(data_list[2])) + ',' + str(float(data_list[3])) + ',' + \
                                  str(float(data_list[4])) + ',' + str(float(data_list[5]))
                        s_insert_data1['data'] = s_data1
                        pca_name = list(s_line_list[0])
                        pca1 = pca_name[0:(len(pca_name) / 2)]
                        pca2 = pca_name[(len(pca_name) / 2):(len(pca_name))]
                        s_insert_data1['x'] = "PC" + str("".join(pca1))
                        s_insert_data1['y'] = "PC" + str("".join(pca2))
                        s_insert_datas1.append(s_insert_data1)
            else:
                # 有分组时
                for s_line in s_data[1:]:
                    s_line_list = s_line.strip().split('\t')
                    for n in range(1, len(s_line_list0)):
                        s_insert_data1 = {
                            "pca_id": self.check_objectid(main_id),
                            "type": "ellipse",
                            "name": s_line_list0[n],
                            "group": s_line_list0[n],
                            "method": "group"
                        }
                        data_list = s_line_list[n].split(",")
                        if len(data_list) == 6:
                            s_data1 = str(float(data_list[0])) + ',' + str(float(data_list[1])) + ',' + \
                                      str(float(data_list[2])) + ',' + str(float(data_list[3])) + ',' + \
                                      str(float(data_list[4])) + ',' + str(float(data_list[5]))
                            s_insert_data1['data'] = s_data1
                            pca_name = list(s_line_list[0])
                            pca1 = pca_name[0:(len(pca_name) / 2)]
                            pca2 = pca_name[(len(pca_name) / 2):(len(pca_name))]
                            s_insert_data1['x'] = "PC" + str("".join(pca1))
                            s_insert_data1['y'] = "PC" + str("".join(pca2))
                            s_insert_datas1.append(s_insert_data1)
            if len(s_insert_datas1) != 0:
                self.col_insert_data("pca_ellipse_detail", s_insert_datas1)
        g_pca_data_insert = {
            "name": "name",
            "data": 'data',
            "group": "group",
            "data_option": ["group", "sample"],
            "condition": {'type': "ellipse"}
        }
        g_update_dict = {
            "ellipse_data": json.dumps(g_pca_data_insert)
        }
        self.update_db_record("pca", {"_id": self.check_objectid(main_id)}, g_update_dict)

    def add_pca_box_detail(self, main_id, output_dir):
        """
        导入画盒形图数据：pca_site.xls
        """
        target_file = os.path.join(output_dir, 'pca_sites.xls')
        if os.path.exists(os.path.join(output_dir, 'group.xls')):
            group_file = os.path.join(output_dir, 'group.xls')
            site = pd.read_table(target_file)
            group = pd.read_table(group_file)
            site.rename(columns={site.columns[0]: "name"}, inplace=True)
            group.rename(columns={group.columns[0]: "name", group.columns[1]: "group"}, inplace=True)
            g_key = group.columns[0]
            new_site = pd.merge(site, group, on=g_key)
            new_site.to_csv(os.path.join(output_dir, 'pca_sites_group.xls'), sep='\t', index=0)
            new_target_file = os.path.join(output_dir, 'pca_sites_group.xls')
            box_df = pd.read_table(new_target_file)
            col_name = box_df.columns
            c = col_name[1]
            box_df1 = box_df[['name', c, 'group']]
            box_df2 = box_df1.groupby('group').describe()
            column_dict = {}
            for n in box_df2[c].columns:
                if 'max' in n:
                    column_dict['Maximum'] = n
                elif 'min' in n:
                    column_dict['Minimum'] = n
                elif '25%' in n:
                    column_dict['Upper quartile'] = n
                elif '75%' in n:
                    column_dict['Lower quartile'] = n
            new_box_df2 = box_df2[c][column_dict.values()]
            new_box_df2.columns = column_dict.keys()
            box_df3 = box_df1.groupby('group').median()
            new_box_df3 = box_df3[c].to_frame()
            new_box_df3.columns = ['Median']
            new_box_df = pd.concat([new_box_df2, new_box_df3], axis=1)
            new_box_df.columns = ['min', 'q3', 'q1', 'max', "median"]
            new_box_df.to_csv(os.path.join(output_dir, c+'_box_data.xls'), sep='\t')
            new_target_file = os.path.join(output_dir, c+'_box_data.xls')
            with open(new_target_file, 'r') as r:
                data = r.readlines()
                pca_data = data[0].strip().split('\t')
                insert_datas = []
                for line in data[1:]:
                    line_list = line.strip().split('\t')
                    insert_data = {
                        "pca_id": self.check_objectid(main_id),
                        "type": "box"
                    }
                    for j in range(len(pca_data)):
                        try:
                            float(line_list[j])
                        except:
                            s_data = line_list[j]
                        else:
                            if float(line_list[j]) < 0.000001:
                                s_data = '{:e}'.format(float(line_list[j]))
                            else:
                                s_data = float(line_list[j])
                        insert_data[pca_data[j]] = s_data
                    insert_data["name"] = insert_data["group"]
                    insert_datas.append(insert_data)
            self.col_insert_data("pca_box_detail", insert_datas)
        else:
            target_file1 = target_file
            box_df = pd.read_table(target_file1)
            box_df1 = box_df.describe().T
            column_dict = {}
            for n in box_df1.columns:
                if 'max' in n:
                    column_dict['Maximum'] = n
                elif 'min' in n:
                    column_dict['Minimum'] = n
                elif '25%' in n:
                    column_dict['Upper quartile'] = n
                elif '75%' in n:
                    column_dict['Lower quartile'] = n
            box_df2 = box_df1[column_dict.values()]
            box_df2.columns = column_dict.keys()
            box_df3 = box_df.median().to_frame()
            box_df4 = pd.concat([box_df2, box_df3], axis=1)
            box_df4.columns = ['min', 'q3', 'q1', 'max', "median"]
            box_df4.to_csv(os.path.join(output_dir, 'all_box_data.xls'), sep='\t')
            new_target_file = os.path.join(output_dir, 'all_box_data.xls')
            with open(new_target_file, 'r') as r:
                insert_datas = []
                data = r.readlines()
                pca_data = data[0].strip().split('\t')
                insert_data = {
                    "pca_id": self.check_objectid(main_id),
                    "type": "box",
                    "name": "all"
                }
                line = data[1]
                line_list = line.strip().split('\t')
                for j in range(len(pca_data)):
                    try:
                        float(line_list[j+1])
                    except:
                        s_data = line_list[j+1]
                    else:
                        if float(line_list[j+1]) < 0.000001:
                            s_data = '{:e}'.format(float(line_list[j+1]))
                        else:
                            s_data = float(line_list[j+1])
                    insert_data[pca_data[j]] = s_data
                insert_data['group'] = 'all'
                insert_datas.append(insert_data)
            self.col_insert_data("pca_box_detail", insert_datas)
        pca_data_insert = {
            "name": "name",
            "category": "group",
            "condition": {'type': "box"}
        }
        update_dict = {
            "box_data": json.dumps(pca_data_insert)
        }
        self.update_db_record("pca", {"_id": self.check_objectid(main_id)}, update_dict)

    def add_pca_result_detail(self, main_id, output_dir):
        """
        导入pca表格数据:导入PCA分析结果表
        """
        target_file1 = os.path.join(output_dir, 'pca_sites.xls')
        with open(target_file1, 'r') as s:
            pca_data_insert = []
            data1 = s.readlines()
            pca_data1 = data1[0].strip().split('\t')
            for line1 in data1[1:]:
                line_list1 = line1.strip().split('\t')
                insert_datas1 = []
                table_dict = {
                    "pca_id": self.check_objectid(main_id),
                    "type": "table"
                }
                for j in range(len(pca_data1)):
                    try:
                        float(line_list1[j])
                    except:
                        s_data1 = line_list1[j]
                    else:
                        if float(line_list1[j]) < 0.000001:
                            s_data1 = '{:e}'.format(float(line_list1[j]))
                        else:
                            s_data1 = float(line_list1[j])
                    table_dict[pca_data1[j]] = s_data1
                insert_datas1.append(table_dict)
                self.col_insert_data("pca_result_detail", insert_datas1)
            for j in range(len(pca_data1)):
                if j == 0:
                    pca_insert_dict = {
                        "field": pca_data1[j],
                        "title": "",
                        "sort": "false",
                        "filter": "false",
                        "type": "string",
                    }
                else:
                    pca_insert_dict = {
                        "field": pca_data1[j],
                        "title": pca_data1[j],
                        "sort": "false",
                        "filter": "false",
                        "type": "float",
                    }
                pca_data_insert.append(pca_insert_dict)
        pca_dict_insert = {
            "column": pca_data_insert,
            "condition": {'type': "table"}
        }
        update_dict = {
            "PCA_result_table": json.dumps(pca_dict_insert)
        }
        self.update_db_record("pca", {"_id": self.check_objectid(main_id)}, update_dict)

    def add_pca_importance_detail(self, main_id, output_dir):
        """
        导入pca表格数据:导入主成分解释表
        """
        target_file2 = os.path.join(output_dir, 'pca_importance.xls')
        axis_data = []
        with open(target_file2, 'r') as s:
            pca_data_insert = []
            data1 = s.readlines()
            pca_data1 = data1[0].strip().split('\t')
            for line1 in data1[1:]:
                line_list1 = line1.strip().split('\t')
                insert_datas1 = []
                table_dict = {
                    "pca_id": self.check_objectid(main_id),
                    "type": "table"
                }
                for j in range(len(pca_data1)):
                    try:
                        float(line_list1[j])
                    except:
                        s_data1 = line_list1[j]
                    else:
                        s_data1 = float(line_list1[j])
                    table_dict[pca_data1[j]] = s_data1
                axis_data.append(line_list1[0] + "(" + str(float(line_list1[1])*100) + "%)")
                insert_datas1.append(table_dict)
                self.col_insert_data("pca_interpretation_detail", insert_datas1)
            for j in range(len(pca_data1)):
                if j == 0:
                    pca_insert_dict = {
                        "field": pca_data1[j],
                        "title": "",
                        "sort": "false",
                        "filter": "false",
                        "type": "string",
                    }
                else:
                    pca_insert_dict = {
                        "field": pca_data1[j],
                        "title": pca_data1[j],
                        "sort": "false",
                        "filter": "false",
                        "type": "float",
                    }
                pca_data_insert.append(pca_insert_dict)
        pca_dict_insert = {
            "column": pca_data_insert,
            "condition": {'type': "table"}
        }
        axis_data_insert = {
            "data": axis_data
        }
        update_dict = {
            "PCA_interpretation_table": json.dumps(pca_dict_insert),
            "axis_data": json.dumps(axis_data_insert)
        }
        self.update_db_record("pca", {"_id": self.check_objectid(main_id)}, update_dict)

    def add_pca_rotation_detail(self, main_id, output_dir):
        """
        导入pca表格数据:导入主成分贡献度表
        """
        target_file3 = os.path.join(output_dir, 'pca_rotation_all.xls')
        with open(target_file3, 'r') as s:
            pca_data_insert = []
            data1 = s.readlines()
            pca_data1 = data1[0].strip().split('\t')
            for line1 in data1[1:]:
                line_list1 = line1.strip().split('\t')
                insert_datas1 = []
                table_dict = {
                    "pca_id": self.check_objectid(main_id),
                    "type": "table"
                }
                for j in range(len(pca_data1)):
                    try:
                        float(line_list1[j])
                    except:
                        s_data1 = line_list1[j]
                    else:
                        s_data1 = float(line_list1[j])
                    table_dict[pca_data1[j]] = s_data1
                insert_datas1.append(table_dict)
                self.col_insert_data("pca_contribution_detail", insert_datas1)
            for j in range(len(pca_data1)):
                if j == 0:
                    pca_insert_dict = {
                        "field": pca_data1[j],
                        "title": "Sample_ID",
                        "sort": "false",
                        "filter": "false",
                        "type": "string",
                    }
                else:
                    pca_insert_dict = {
                        "field": pca_data1[j],
                        "title": pca_data1[j],
                        "sort": "false",
                        "filter": "false",
                        "type": "float",
                    }
                pca_data_insert.append(pca_insert_dict)
        pca_dict_insert = {
            "column": pca_data_insert,
            "condition": {'type': "table"}
        }
        update_dict = {
            "PCA_contribution_table": json.dumps(pca_dict_insert)
        }
        self.update_db_record("pca", {"_id": self.check_objectid(main_id)}, update_dict)


if __name__ == "__main__":
    a = Pca(None)
    a.add_pca_detail("5ec3490000b2bf4f5912d5ee", "/mnt/ilustre/users/sanger-dev/workspace/20200901/Pca_pca20200901090135/output/pca")
    a.add_pca_ellipse_detail("5ec3490000b2bf4f5912d5ee", "/mnt/ilustre/users/sanger-dev/workspace/20200901/Pca_pca20200901090135/output/pca")
    a.add_pca_box_detail("5ec3490000b2bf4f5912d5ee", "/mnt/ilustre/users/sanger-dev/workspace/20200901/Pca_pca20200901090135/output/pca")
    a.add_pca_result_detail("5ec3490000b2bf4f5912d5ee", "/mnt/ilustre/users/sanger-dev/workspace/20200901/Pca_pca20200901090135/output/pca")
    a.add_pca_importance_detail("5ec3490000b2bf4f5912d5ee", "/mnt/ilustre/users/sanger-dev/workspace/20200901/Pca_pca20200901090135/output/pca")
    a.add_pca_rotation_detail("5ec3490000b2bf4f5912d5ee", "/mnt/ilustre/users/sanger-dev/workspace/20200901/Pca_pca20200901090135/output/pca")
