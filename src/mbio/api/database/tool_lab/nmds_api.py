# !usr/bin/python
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


class NmdsApi(ApiBase):
    def __init__(self, bind_object):
        super(NmdsApi, self).__init__(bind_object)
        self._db_name = 'sanger_tool_lab'

    def add_nmds_detail(self, main_id, output_dir):
        """
        导入nmds画散点图数据
        """
        target_file0 = os.path.join(output_dir, 'nmds_sites.xls')
        re_site = pd.read_table(target_file0)
        re_site.rename(columns={re_site.columns[0]: "name", "MDS1": "NMDS1", "MDS2": "NMDS2"}, inplace=True)
        re_site.to_csv(os.path.join(output_dir, 'new_nmds_sites.xls'), sep='\t', index=0)
        target_file = os.path.join(output_dir, 'new_nmds_sites.xls')
        if os.path.exists(os.path.join(output_dir, 'group.xls')):
            group_file = os.path.join(output_dir, 'group.xls')
            site = pd.read_table(target_file)
            group = pd.read_table(group_file)
            group.rename(columns={group.columns[0]: "name", group.columns[1]: "group"}, inplace=True)
            g_key = group.columns[0]
            new_site = pd.merge(site, group, on=g_key)
            new_site.to_csv(os.path.join(output_dir, 'nmds_sites_group.xls'), sep='\t', index=0)
            new_target_file = os.path.join(output_dir, 'nmds_sites_group.xls')
        else:
            new_target_file = target_file
        with open(new_target_file, 'r') as r:
            insert_datas = []
            data = r.readlines()
            nmds_data = data[0].strip().split('\t')
            for line in data[1:]:
                line_list = line.strip().split('\t')
                insert_data = {
                    "nmds_id": self.check_objectid(main_id),
                    "type": "scatter",
                }
                for j in range(len(nmds_data)):
                    try:
                        float(line_list[j])
                    except:
                        s_data = line_list[j]
                    else:
                        s_data = float(line_list[j])
                    insert_data[nmds_data[j]] = s_data
                    if 'group' not in nmds_data:
                        insert_data['group'] = ''
                insert_datas.append(insert_data)
        self.col_insert_data("nmds_detail", insert_datas)
        new_nmds_data = nmds_data
        if 'name' in new_nmds_data:
            new_nmds_data.remove('name')
        if 'group' in new_nmds_data:
            new_nmds_data.remove('group')
        nmds_data_insert = {
            "name": "name",
            "data": new_nmds_data,
            "category": "group",
            "condition": {'type': "scatter"}
        }
        update_dict = {
            "scatter_data": json.dumps(nmds_data_insert),
        }
        self.update_db_record("nmds", {"_id": self.check_objectid(main_id)}, update_dict)

    def add_nmds_ellipse_detail(self, main_id, output_dir):
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
                        "nmds_id": self.check_objectid(main_id),
                        "type": "ellipse",
                        "name": s_line_list0[1],
                        "group": s_line_list0[1],
                        "method": "sample"
                    }
                    data_list = s_line_list[1].split(",")
                    if len(data_list) == 6:
                        # 数据出来为 m1, c11, c12, m2, c21, c22
                        # s_data1 = str(float(data_list[5])) + ',' + str(float(data_list[4])) + ',' + \
                        #           str(float(data_list[3])) + ',' + str(float(data_list[0])) + ',' + \
                        #           str(float(data_list[2])) + ',' + str(float(data_list[1]))
                        s_data1 = str(float(data_list[0])) + ',' + str(float(data_list[1])) + ',' + \
                                  str(float(data_list[2])) + ',' + str(float(data_list[3])) + ',' + \
                                  str(float(data_list[4])) + ',' + str(float(data_list[5]))
                        s_insert_data1['data'] = s_data1
                        nmds_name = list(s_line_list[0])
                        ndms1 = nmds_name[0:(len(nmds_name)/2)]
                        nmds2 = nmds_name[(len(nmds_name)/2):(len(nmds_name))]
                        s_insert_data1['x'] = "NMDS" + str("".join(ndms1))
                        s_insert_data1['y'] = "NMDS" + str("".join(nmds2))
                        s_insert_datas1.append(s_insert_data1)
            else:
                # 有分组时
                for s_line in s_data[1:]:
                    s_line_list = s_line.strip().split('\t')
                    for n in range(1, len(s_line_list0)):
                        s_insert_data1 = {
                            "nmds_id": self.check_objectid(main_id),
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
                            nmds_name = list(s_line_list[0])
                            ndms1 = nmds_name[0:(len(nmds_name) / 2)]
                            nmds2 = nmds_name[(len(nmds_name) / 2):(len(nmds_name))]
                            s_insert_data1['x'] = "NMDS" + str("".join(ndms1))
                            s_insert_data1['y'] = "NMDS" + str("".join(nmds2))
                            s_insert_datas1.append(s_insert_data1)
            if len(s_insert_datas1) != 0:
                self.col_insert_data("nmds_ellipse_detail", s_insert_datas1)
        g_nmds_data_insert = {
            "name": "name",
            "data": 'data',
            "group": "group",
            "data_option": ["group", "sample"],
            "condition": {'type': "ellipse"}
        }
        g_update_dict = {
            "ellipse_data": json.dumps(g_nmds_data_insert)
        }
        self.update_db_record("nmds", {"_id": self.check_objectid(main_id)}, g_update_dict)

    def add_nmds_table_detail(self, main_id, output_dir):
        """
        导入nmds表格数据,nmds_sites.xls
        """
        target_file = os.path.join(output_dir, 'new_nmds_sites.xls')
        with open(target_file, 'r') as s:
            nmds_data_insert = []
            data1 = s.readlines()
            nmds_data1 = data1[0].strip().split('\t')
            for line1 in data1[1:]:
                line_list1 = line1.strip().split('\t')
                insert_datas1 = []
                table_dict = {
                    "nmds_id": self.check_objectid(main_id),
                    "type": "table"
                }
                for j in range(len(nmds_data1)):
                    try:
                        float(line_list1[j])
                    except:
                        s_data1 = line_list1[j]
                    else:
                        s_data1 = float(line_list1[j])
                    table_dict[nmds_data1[j]] = s_data1
                insert_datas1.append(table_dict)
                self.col_insert_data("nmds_result_detail", insert_datas1)
            for j in range(len(nmds_data1)):
                if j == 0:
                    nmds_insert_dict = {
                        "field": nmds_data1[j],
                        "title": "",
                        "sort": "false",
                        "filter": "false",
                        "type": "string",
                    }
                else:
                    nmds_insert_dict = {
                        "field": nmds_data1[j],
                        "title": nmds_data1[j],
                        "sort": "false",
                        "filter": "false",
                        "type": "float",
                    }
                nmds_data_insert.append(nmds_insert_dict)
        nmds_dict_insert = {
            "column": nmds_data_insert,
            "condition": {'type': "table"}
        }
        update_dict = {
            "nmds_result_table": json.dumps(nmds_dict_insert)
        }
        self.update_db_record("nmds", {"_id": self.check_objectid(main_id)}, update_dict)

    def add_nmds_stress_detail(self, main_id, output_dir):
        """
        导入nmds画stress数据
        """
        stress_file = os.path.join(output_dir, 'nmds_stress.xls')
        insert_stress_datas = []
        with open(stress_file, 'rb') as s:
            i = 0
            for line in s:
                if i == 0:
                    i = 1
                else:
                    stress = line.strip('\n')
        insert_stress_data = {
            "nmds_id": self.check_objectid(main_id),
            "type": "text",
            'text': 'stress:'+stress,
            'name': 'stress'
        }
        insert_stress_datas.append(insert_stress_data)
        self.col_insert_data("nmds_detail", insert_stress_datas)
        nmds_data_insert = {
            "name": "name",
            "data": 'text',
            "condition": {'type': "text"}
        }
        update_dict = {
            "text_data": json.dumps(nmds_data_insert)
        }
        self.update_db_record("nmds", {"_id": self.check_objectid(main_id)}, update_dict)


if __name__ == "__main__":
    a = NmdsApi(None)
    # a.add_nmds_detail("5e9eb43b17b2bf4b84491222", "/mnt/ilustre/users/sanger-dev/workspace/20200713/Nmds_a25dprj7o13r28rjp8uhnq46ld_0713103447913660_7217/output/nmds")
    a.add_nmds_ellipse_detail("5e9eb43b17b2bf4b81111222", "/mnt/ilustre/users/sanger-dev/workspace/20200713/Nmds_a25dprj7o13r28rjp8uhnq46ld_0713103447913660_7217/output/nmds")
    # a.add_nmds_table_detail("5e9eb43b17b2bf4b84491222", "/mnt/ilustre/users/sanger-dev/workspace/20200713/Nmds_a25dprj7o13r28rjp8uhnq46ld_0713103447913660_7217/output/nmds")
    # a.add_nmds_stress_detail("5e9eb43b17b2bf4b84491222", "/mnt/ilustre/users/sanger-dev/workspace/20200713/Nmds_a25dprj7o13r28rjp8uhnq46ld_0713103447913660_7217/output/nmds")
