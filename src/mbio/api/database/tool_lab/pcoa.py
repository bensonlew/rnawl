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


class Pcoa(ApiBase):
    def __init__(self, bind_object):
        super(Pcoa, self).__init__(bind_object)
        self._db_name = 'sanger_tool_lab'

    def add_pcoa_detail(self, main_id, output_dir):
        """
        导入pcoa画散点图数据
        """
        target_file = os.path.join(output_dir, 'pcoa_sites.xls')
        if os.path.exists(os.path.join(output_dir, 'group.xls')):
            group_file = os.path.join(output_dir, 'group.xls')
            site = pd.read_table(target_file)
            group = pd.read_table(group_file)
            site.rename(columns={site.columns[0]: "name"}, inplace=True)
            group.rename(columns={group.columns[0]: "name", group.columns[1]: "group"}, inplace=True)
            g_key = group.columns[0]
            new_site = pd.merge(site, group, on=g_key)
            new_site.to_csv(os.path.join(output_dir, 'pcoa_sites_group.xls'), sep='\t', index=0)
            new_target_file = os.path.join(output_dir, 'pcoa_sites_group.xls')
        else:
            site = pd.read_table(target_file)
            site.rename(columns={site.columns[0]: "name"}, inplace=True)
            site.to_csv(os.path.join(output_dir, 'pcoa_sites_new.xls'), sep='\t', index=0)
            new_target_file = os.path.join(output_dir, 'pcoa_sites_new.xls')
        with open(new_target_file, 'r') as r:
            insert_datas = []
            data = r.readlines()
            pcoa_data = data[0].strip().split('\t')
            for line in data[1:]:
                line_list = line.strip().split('\t')
                insert_data = {
                    "pcoa_id": self.check_objectid(main_id),
                    "type": "scatter",
                }
                for j in range(len(pcoa_data)):
                    try:
                        float(line_list[j])
                    except:
                        s_data = line_list[j]
                    else:
                        if float(line_list[j]) < 0.000001:
                            s_data = '{:e}'.format(float(line_list[j]))
                        else:
                            s_data = float(line_list[j])
                    insert_data[pcoa_data[j]] = s_data
                    if 'group' not in pcoa_data:
                        insert_data['group'] = ''
                insert_datas.append(insert_data)
        self.col_insert_data("pcoa_detail", insert_datas)
        new_pcoa_data = pcoa_data
        if 'name' in new_pcoa_data:
            new_pcoa_data.remove('name')
        if 'group' in new_pcoa_data:
            new_pcoa_data.remove('group')
        pcoa_data_insert = {
            "name": "name",
            "data": new_pcoa_data,
            "category": "group",
            "condition": {'type': "scatter"}
        }
        update_dict = {
            "scatter_data": json.dumps(pcoa_data_insert),
        }
        self.update_db_record("pcoa", {"_id": self.check_objectid(main_id)}, update_dict)

    def add_pcoa_ellipse_detail(self, main_id, output_dir):
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
                        "pcoa_id": self.check_objectid(main_id),
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
                        pcoa_name = list(s_line_list[0])
                        pcoa1 = pcoa_name[0:(len(pcoa_name) / 2)]
                        pcoa2 = pcoa_name[(len(pcoa_name) / 2):(len(pcoa_name))]
                        s_insert_data1['x'] = "PC" + str("".join(pcoa1))
                        s_insert_data1['y'] = "PC" + str("".join(pcoa2))
                        s_insert_datas1.append(s_insert_data1)
            else:
                # 有分组时
                for s_line in s_data[1:]:
                    s_line_list = s_line.strip().split('\t')
                    for n in range(1, len(s_line_list0)):
                        s_insert_data1 = {
                            "pcoa_id": self.check_objectid(main_id),
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
                            pcoa_name = list(s_line_list[0])
                            pcoa1 = pcoa_name[0:(len(pcoa_name) / 2)]
                            pcoa2 = pcoa_name[(len(pcoa_name) / 2):(len(pcoa_name))]
                            s_insert_data1['x'] = "PC" + str("".join(pcoa1))
                            s_insert_data1['y'] = "PC" + str("".join(pcoa2))
                            s_insert_datas1.append(s_insert_data1)
            if len(s_insert_datas1) != 0:
                self.col_insert_data("pcoa_ellipse_detail", s_insert_datas1)
        g_pcoa_data_insert = {
            "name": "name",
            "data": 'data',
            "group": "group",
            "data_option": ["group", "sample"],
            "condition": {'type': "ellipse"}
        }
        g_update_dict = {
            "ellipse_data": json.dumps(g_pcoa_data_insert)
        }
        self.update_db_record("pcoa", {"_id": self.check_objectid(main_id)}, g_update_dict)

    def add_pcoa_box_detail(self, main_id, output_dir):
        """
        导入画盒形图数据：pcoa_sites_group.xls
        """
        target_file = os.path.join(output_dir, 'pcoa_sites.xls')
        if os.path.exists(os.path.join(output_dir, 'group.xls')):
            group_file = os.path.join(output_dir, 'group.xls')
            site = pd.read_table(target_file)
            group = pd.read_table(group_file)
            site.rename(columns={site.columns[0]: "name"}, inplace=True)
            group.rename(columns={group.columns[0]: "name", group.columns[1]: "group"}, inplace=True)
            g_key = group.columns[0]
            new_site = pd.merge(site, group, on=g_key)
            new_site.to_csv(os.path.join(output_dir, 'pcoa_sites_group.xls'), sep='\t', index=0)
            new_target_file = os.path.join(output_dir, 'pcoa_sites_group.xls')
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
            new_target_file = os.path.join(output_dir, c + '_box_data.xls')
            with open(new_target_file, 'r') as r:
                data = r.readlines()
                insert_datas = []
                pcoa_data = data[0].strip().split('\t')
                for line in data[1:]:
                    insert_data = {
                        "pcoa_id": self.check_objectid(main_id),
                        "type": "box"
                    }
                    line_list = line.strip().split('\t')
                    for j in range(len(pcoa_data)):
                        try:
                            float(line_list[j])
                        except:
                            s_data = line_list[j]
                        else:
                            if float(line_list[j]) < 0.000001:
                                s_data = '{:e}'.format(float(line_list[j]))
                            else:
                                s_data = float(line_list[j])
                        insert_data[pcoa_data[j]] = s_data
                    insert_data["name"] = insert_data["group"]
                    insert_datas.append(insert_data)
            self.col_insert_data("pcoa_box_detail", insert_datas)
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
                pcoa_data = data[0].strip().split('\t')
                insert_data = {
                    "pcoa_id": self.check_objectid(main_id),
                    "type": "box",
                    "name": "all"
                }
                line = data[1]
                line_list = line.strip().split('\t')
                for j in range(len(pcoa_data)):
                    try:
                        float(line_list[j+1])
                    except:
                        s_data = line_list[j+1]
                    else:
                        if float(line_list[j+1]) < 0.000001:
                            s_data = '{:e}'.format(float(line_list[j+1]))
                        else:
                            s_data = float(line_list[j+1])
                    insert_data[pcoa_data[j]] = s_data
                insert_data['group'] = 'all'
                insert_datas.append(insert_data)
            self.col_insert_data("pcoa_box_detail", insert_datas)
        pcoa_data_insert = {
            "name": "name",
            "category": "group",
            "condition": {'type': "box"}
        }
        update_dict = {
            "box_data": json.dumps(pcoa_data_insert)
        }
        self.update_db_record("pcoa", {"_id": self.check_objectid(main_id)}, update_dict)

    def add_pcoa_result_detail(self, main_id, output_dir):
        """
        导入pcoa表格数据:导入PCoA分析结果表, pcoa_sites.xls
        """
        target_file1 = os.path.join(output_dir, 'pcoa_sites.xls')
        with open(target_file1, 'r') as s:
            pcoa_data_insert = []
            data1 = s.readlines()
            pcoa_data1 = data1[0].strip().split('\t')
            for line1 in data1[1:]:
                line_list1 = line1.strip().split('\t')
                insert_datas1 = []
                table_dict = {
                    "pcoa_id": self.check_objectid(main_id),
                    "type": "table"
                }
                for j in range(len(pcoa_data1)):
                    try:
                        float(line_list1[j])
                    except:
                        s_data1 = line_list1[j]
                    else:
                        if float(line_list1[j]) < 0.000001:
                            s_data1 = '{:e}'.format(float(line_list1[j]))
                        else:
                            s_data1 = float(line_list1[j])
                    table_dict[pcoa_data1[j]] = s_data1
                insert_datas1.append(table_dict)
                self.col_insert_data("pcoa_result_detail", insert_datas1)
            for j in range(len(pcoa_data1)):
                if j == 0:
                    pcoa_insert_dict = {
                        "field": pcoa_data1[j],
                        "title": "",
                        "sort": "false",
                        "filter": "false",
                        "type": "string",
                    }
                else:
                    pcoa_insert_dict = {
                        "field": pcoa_data1[j],
                        "title": pcoa_data1[j],
                        "sort": "false",
                        "filter": "false",
                        "type": "float",
                    }
                pcoa_data_insert.append(pcoa_insert_dict)
        pcoa_dict_insert = {
            "column": pcoa_data_insert,
            "condition": {'type': "table"}
        }
        update_dict = {
            "PCoA_result_table": json.dumps(pcoa_dict_insert)
        }
        self.update_db_record("pcoa", {"_id": self.check_objectid(main_id)}, update_dict)

    def add_pcoa_eigenvalues_detail(self, main_id, output_dir):
        """
        导入pcoa表格数据:导入矩阵特征值表, pcoa_eigenvalues.xls
        """
        target_file2 = os.path.join(output_dir, 'pcoa_eigenvalues.xls')
        with open(target_file2, 'r') as s:
            pcoa_data_insert = []
            data1 = s.readlines()
            pcoa_data0 = data1[0].strip().split('\t')
            pcoa_data1 = ['name', pcoa_data0[0]]
            for line1 in data1[1:]:
                line_list1 = line1.strip().split('\t')
                insert_datas1 = []
                table_dict = {
                    "pcoa_id": self.check_objectid(main_id),
                    "type": "table"
                }
                for j in range(len(pcoa_data1)):
                    try:
                        float(line_list1[j])
                    except:
                        s_data1 = line_list1[j]
                    else:
                        s_data1 = float(line_list1[j])
                    table_dict[pcoa_data1[j]] = s_data1
                insert_datas1.append(table_dict)
                self.col_insert_data("pcoa_eigenvalues_detail", insert_datas1)
            for j in range(len(pcoa_data1)):
                if j == 0:
                    pcoa_insert_dict = {
                        "field": pcoa_data1[j],
                        "title": "",
                        "sort": "false",
                        "filter": "false",
                        "type": "string",
                    }
                else:
                    pcoa_insert_dict = {
                        "field": pcoa_data1[j],
                        "title": 'eigenvalues',
                        "sort": "false",
                        "filter": "false",
                        "type": "float",
                    }
                pcoa_data_insert.append(pcoa_insert_dict)
        pcoa_dict_insert = {
            "column": pcoa_data_insert,
            "condition": {'type': "table"}
        }
        update_dict = {
            "PCoA_eigenvalue_table": json.dumps(pcoa_dict_insert)
        }
        self.update_db_record("pcoa", {"_id": self.check_objectid(main_id)}, update_dict)

    def add_pcoa_eigenvaluespre_detail(self, main_id, output_dir):
        """
        导入pcoa表格数据:导入主成分解释表, pcoa_eigenvaluespre.xls
        """
        target_file3 = os.path.join(output_dir, 'pcoa_eigenvaluespre.xls')
        axis_data = []
        with open(target_file3, 'r') as s:
            pcoa_data_insert = []
            data1 = s.readlines()
            pcoa_data0 = data1[0].strip().split('\t')
            pcoa_data1 = ['name', pcoa_data0[0]]
            for line1 in data1[1:]:
                line_list1 = line1.strip().split('\t')
                insert_datas1 = []
                table_dict = {
                    "pcoa_id": self.check_objectid(main_id),
                    "type": "table"
                }
                for j in range(len(pcoa_data1)):
                    try:
                        float(line_list1[j])
                    except:
                        s_data1 = line_list1[j]
                    else:
                        s_data1 = float(line_list1[j])
                    table_dict[pcoa_data1[j]] = s_data1
                axis_data.append(line_list1[0] + "(" + str(float(line_list1[1])*100) + "%)")
                insert_datas1.append(table_dict)
                self.col_insert_data("pcoa_interpretation_detail", insert_datas1)
            for j in range(len(pcoa_data1)):
                if j == 0:
                    pcoa_insert_dict = {
                        "field": pcoa_data1[j],
                        "title": "",
                        "sort": "false",
                        "filter": "false",
                        "type": "string",
                    }
                else:
                    pcoa_insert_dict = {
                        "field": pcoa_data1[j],
                        "title": "Proportion of Variance",
                        "sort": "false",
                        "filter": "false",
                        "type": "float",
                    }
                pcoa_data_insert.append(pcoa_insert_dict)
        pcoa_dict_insert = {
            "column": pcoa_data_insert,
            "condition": {'type': "table"}
        }
        axis_data_insert = {
            "data": axis_data
        }
        update_dict = {
            "PCoA_interpretation_table": json.dumps(pcoa_dict_insert),
            "axis_data": json.dumps(axis_data_insert)
        }
        self.update_db_record("pcoa", {"_id": self.check_objectid(main_id)}, update_dict)

    def add_pcoa_distance_detail(self, main_id, output_dir):
        """
        导入pcoa表格数据:导入PCoA样本间距离矩阵表, distance.xls
        """
        target_file1 = os.path.join(output_dir, 'distance.xls')
        with open(target_file1, 'r') as s:
            pcoa_data_insert = []
            data1 = s.readlines()
            pcoa_data1 = data1[0].strip().split('\t')
            pcoa_data1.insert(0, 'sample')
            for line1 in data1[1:]:
                line_list1 = line1.strip().split('\t')
                insert_datas1 = []
                table_dict = {
                    "pcoa_id": self.check_objectid(main_id),
                    "type": "table"
                }
                table_dict['name'] = pcoa_data1[0]
                for j in range(len(pcoa_data1)):
                    try:
                        float(line_list1[j])
                    except:
                        s_data1 = line_list1[j]
                    else:
                        if float(line_list1[j]) < 0.000001:
                            s_data1 = '{:e}'.format(float(line_list1[j]))
                        else:
                            s_data1 = float(line_list1[j])
                    table_dict[pcoa_data1[j]] = s_data1
                insert_datas1.append(table_dict)
                self.col_insert_data("pcoa_distance_detail", insert_datas1)
            for j in range(len(pcoa_data1)):
                if j == 0:
                    pcoa_insert_dict = {
                        "field": pcoa_data1[j],
                        "title": "",
                        "sort": "false",
                        "filter": "false",
                        "type": "string",
                    }
                else:
                    pcoa_insert_dict = {
                        "field": pcoa_data1[j],
                        "title": pcoa_data1[j],
                        "sort": "false",
                        "filter": "false",
                        "type": "float",
                    }
                pcoa_data_insert.append(pcoa_insert_dict)
        pcoa_dict_insert = {
            "column": pcoa_data_insert,
            "condition": {'type': "table"}
        }
        update_dict = {
            "PCoA_distance_table": json.dumps(pcoa_dict_insert)
        }
        self.update_db_record("pcoa", {"_id": self.check_objectid(main_id)}, update_dict)

if __name__ == "__main__":
    a = Pcoa(None)
    a.add_pcoa_detail("5ea91cbf17b2bf65274c3333", "/mnt/ilustre/users/sanger-dev/i-sanger_workspace2/20210316/Pcoa_8u02_86l8u01o50apgpd353mhf0_0316132435807469_9324/output/pcoa")
    a.add_pcoa_ellipse_detail("5ea91cbf17b2bf65274c3333", "/mnt/ilustre/users/sanger-dev/i-sanger_workspace2/20210316/Pcoa_8u02_86l8u01o50apgpd353mhf0_0316132435807469_9324/output/pcoa")
    a.add_pcoa_box_detail("5ea91cbf17b2bf65274c3333", "/mnt/ilustre/users/sanger-dev/i-sanger_workspace2/20210316/Pcoa_8u02_86l8u01o50apgpd353mhf0_0316132435807469_9324/output/pcoa")
    a.add_pcoa_result_detail("5ea91cbf17b2bf65274c3333", "/mnt/ilustre/users/sanger-dev/i-sanger_workspace2/20210316/Pcoa_8u02_86l8u01o50apgpd353mhf0_0316132435807469_9324/output/pcoa")
    a.add_pcoa_eigenvalues_detail("5ea91cbf17b2bf65274c3333", "/mnt/ilustre/users/sanger-dev/i-sanger_workspace2/20210316/Pcoa_8u02_86l8u01o50apgpd353mhf0_0316132435807469_9324/output/pcoa")
    a.add_pcoa_eigenvaluespre_detail("5ea91cbf17b2bf65274c3333", "/mnt/ilustre/users/sanger-dev/i-sanger_workspace2/20210316/Pcoa_8u02_86l8u01o50apgpd353mhf0_0316132435807469_9324/output/pcoa")
    a.add_pcoa_distance_detail("5ea91cbf17b2bf65274c3333", "/mnt/ilustre/users/sanger-dev/i-sanger_workspace2/20210316/Pcoa_8u02_86l8u01o50apgpd353mhf0_0316132435807469_9324/output/pcoa")
