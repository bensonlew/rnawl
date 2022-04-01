# -*- coding: utf-8 -*-
# __author__ = 'linna'

import json
import datetime
from api_base import ApiBase
from bson.objectid import ObjectId
import os
import pandas as pd


class Regression(ApiBase):
    def __init__(self, bind_object):
        super(Regression, self).__init__(bind_object)

    def add_regression(self, main_id, output_dir, input_dir):
        group_dict = {}
        if os.path.exists(os.path.join(output_dir, 'group.xls')):
            group_file = os.path.join(output_dir, 'group.xls')
            group = pd.read_table(group_file)
            with open(group_file)as f:
                lines = f.readlines()
                for line in lines[1:]:
                    group_dict[line.strip().split("\t")[0]] = line.strip().split("\t")[1]
            origin = pd.read_table(input_dir)
            target = pd.merge(origin, group, on='name')
            target.to_csv(os.path.join(output_dir, 'regression_data_group.xls'), sep='\t', index=0)

            target_file = os.path.join(output_dir, 'regression_data_group.xls')
            result = pd.read_table(target_file, header=0)
            sample_dict_list = result.to_dict('records')
        else:
            target_file = os.path.join(input_dir)
            result = pd.read_table(target_file, header=0)
            sample_dict_list = result.to_dict('records')

        site_file = os.path.join(output_dir, 'regression_data.xls')
        result = pd.read_table(site_file)
        site_dict_list = result.to_dict('records')
        rlt_col = result.columns
        for dict1 in site_dict_list:
            if os.path.exists(os.path.join(output_dir, 'group.xls')):
                dict1["name"] = dict1.pop(rlt_col[0])
                dict1["group"] = group_dict[dict1[rlt_col[0]]]
                dict1["type"] = "scatter"
                dict1["y"] = dict1[rlt_col[1]]
                dict1["x"] = dict1[rlt_col[2]]
            else:
                dict1["group"] = "all"
                dict1["type"] = "scatter"
                dict1["y"] = dict1[rlt_col[1]]
                dict1["x"] = dict1[rlt_col[2]]
            dict1["regression_id"] = self.check_objectid(main_id)
        # main_id = ObjectId(main_id)
        # tag_dict = dict(regression_id=main_id)
        # self.col_insert_data('sg_regression_detail', sample_dict_list)
        self.col_insert_data('sg_regression_scatter_detail', site_dict_list)

        target_file = os.path.join(output_dir, 'regression_message.xls')
        with open(target_file, 'rb') as r:
            i = 0
            for line in r:
                if i == 0:
                    i = 1
                else:
                    line = line.strip('\n')
                    line_data = line.split('\t')
                    line_d = ''
                    # if float('%.3f' % float(line_data[6])) > 0:
                    #     line_d = 'y=' + str('%.3f' % float(line_data[5])) + 'x+' + str('%.3f' % float(line_data[6]))
                    # elif float('%.3f' % float(line_data[6])) < 0:
                    #     line_d = 'y=' + str('%.3f' % float(line_data[5])) + 'x' + str('%.3f' % float(line_data[6]))
                    # elif float('%.3f' % float(line_data[6])) == 0.000:
                    #     line_d = 'y=' + str('%.3f' % float(line_data[5])) + 'x'
                    # 产品需求将第一列做X，第二列做Y
                    b = float(line_data[1]) - float(line_data[2])*(1 / float(line_data[5]))
                    if float('%.3f' % b) > 0:
                        line_d = 'y=' + str('%.3f' % (1 / float(line_data[5]))) + 'x-' + str('%.3f' % b)
                    elif float('%.3f' % b) < 0:
                        line_d = 'y=' + str('%.3f' % (1 / float(line_data[5]))) + 'x' + str('%.3f' % b)
                    elif float('%.3f' % b) == 0.000:
                        line_d = 'y=' + str('%.3f' % (1 / float(line_data[5]))) + 'x'
                    R_2 = '%.2f' % float(line_data[0])
                    xmin = float(line_data[2])
                    ymin = float(line_data[1])
                    xmax = float(line_data[4])
                    ymax = float(line_data[3])
                    regression_equation = line_d

        self.col_insert_data('sg_regression_line_detail', [{
            "x1": xmin,
            "y1": ymin,
            "x2": xmax,
            "y2": ymax,
            "type": "divide_line",
            "name": "",
            "group": "",
            "regression_id": self.check_objectid(main_id)
        }])

        self.col_insert_data('sg_regression_text_detail', [{
            "regression_id": self.check_objectid(main_id),
            "text": regression_equation,
            "name": "",
            "group": "",
            "type": "text",
        }])

        self.col_insert_data('sg_regression_text_detail', [{
            "regression_id": self.check_objectid(main_id),
            "text": "R^2 =" + R_2,
            "name": "",
            "group": "",
            "type": "text",
        }])

        update_dict = {
            "scatter_data": json.dumps({
                "name": "name",
                "data": ["x", "y"],
                "category": "group",
                "condition": {'type': "scatter"}
            }),
            "divide_line_data": json.dumps({
                "name": "name",
                "condition": {'type': "divide_line"}
            }),
            "text_data": json.dumps({
                "name": "name",
                "condition": {'type': "text"}
            })
        }
        self.update_db_record("sg_regression", {"_id": self.check_objectid(main_id)}, update_dict)

if __name__ == '__main__':
   a = Regression(None)
   mainid = "5ecf64f217b2bf12a65e6af6"
   output_dir = '/mnt/lustre/users/sanger-dev/wpm2/workspace/20220120/Regression_ilkn_7128mj2rtmm1u6bbhkpf78_0120085307065507_6594/output/regression'
   input = '/mnt/lustre/users/sanger-dev/wpm2/workspace/20220120/Regression_ilkn_7128mj2rtmm1u6bbhkpf78_0120085307065507_6594/remote_input/tooltable/Regression_data1.txt'
   a.add_regression(mainid, output_dir, input)
