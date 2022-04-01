# -*- coding: utf-8 -*-
# __author__ = 'binbinzhao'

import json
import os
from api_base import ApiBase


class RocAnalysis(ApiBase):
    def __init__(self, bind_object):
        super(RocAnalysis, self).__init__(bind_object)

    def add_sg_roc(self, main_id, file_path, subtype):
        data_list = []
        data_list_scatter = []
        data_list_text = []
        files = os.listdir(file_path)
        if subtype:
            subtype_new = "curveCardinal"
        else:
            subtype_new = "curveLinear"
        for file in files:
            with open(file_path + "/" + file) as f:
                lines = f.readlines()
                for i in range(1, len(lines[2].strip().split(" ")[1:])):
                    insert_data = {
                        "roc_id":  self.check_objectid(main_id),
                        "name": "file",
                        "category": file.split(".")[0],
                        "type": "line",
                        "subtype": subtype_new,
                        "x": round(float(lines[2].strip().split(" ")[i]), 2),
                        "y": round(float(lines[3].strip().split(" ")[i]), 2),
                    }
                    if len(files) == 1:
                        insert_data["ymin"] = float(0)
                        insert_data["ymax"] = round(float(lines[3].strip().split(" ")[i]), 2)
                    data_list.append(insert_data)
                listx = lines[2].strip().split(" ")
                listx[1:len(listx)].sort()
                medainx = (float(listx[len(listx)/2]) + float(listx[~(len(listx)/2)])) / 2
                listy = lines[3].strip().split(" ")
                listy[1:len(listy)].sort()
                medainy = (float(listy[len(listy) / 2]) + float(listy[~(len(listy) / 2)])) / 2
                x_scatter = lines[6].strip().split(" ")
                y_scatter = lines[7].strip().split(" ")
                if len(x_scatter) > 1 and len(y_scatter) > 1:
                    insert_data_scatter = {
                        "roc_id": self.check_objectid(main_id),
                        "x": x_scatter[1],
                        "y": y_scatter[1],
                        "name": lines[5].strip().split(" ")[1] + "(" + x_scatter[1]\
                                + "," + y_scatter[1] + ")",
                        "group": "",
                        "type": "scatter"
                    }
                    data_list_scatter.append(insert_data_scatter)
                if len(files) == 1:
                    insert_data_text = {
                        "roc_id": self.check_objectid(main_id),
                        "type": "text",
                        "x": medainx,
                        "y": medainy,
                        "text": "AUC:" + lines[9].strip().split(" ")[1]
                    }
                    data_list_text.append(insert_data_text)
        if len(data_list) == 0:
            self.bind_object.logger.info("{}文件为空！".format(file_path))
        else:
            self.col_insert_data("sg_roc_line", data_list)
            if len(data_list_scatter) > 0:
                self.col_insert_data("sg_roc_scatter", data_list_scatter)
            if len(files) == 1:
                self.col_insert_data("sg_roc_text", data_list_text)

        if len(files) == 1:   # 增加chart_type字段区分单线和多线  add by wuqin at 20200722
            update_dict = {
                "line_data": json.dumps({
                    "name": "name",
                    "condition": {'type': "line"},
                }),
                "scatter_data": json.dumps({
                    "name": "name",
                    "data": ["x", "y"],
                    "category": "group",
                    "condition": {'type': "scatter"}
                }),
                "text_data": json.dumps({
                    "name": "name",
                    "condition": {'type': "text"}
                }),
                "chart_type": json.dumps({
                    "data": "area"
                }),
                "area_data": json.dumps({
                    "name": "name",
                    "condition": {'type': "line"},
                })
            }
        else:
            update_dict = {
                "line_data": json.dumps({
                    "name": "name",
                    "condition": {'type': "line"},
                }),
                "scatter_data": json.dumps({
                    "name": "name",
                    "data": ["x", "y"],
                    "category": "group",
                    "condition": {'type': "scatter"}
                }),
                "text_data": json.dumps({
                    "name": "name",
                    "condition": {'type': "text"}
                }),
                "chart_type": json.dumps({
                    "data": "line"
                }),
                "area_data": json.dumps({
                    "name": "name",
                    "condition": {'type': "line"},
                })
            }
        self.update_db_record("sg_roc", {"_id": self.check_objectid(main_id)}, update_dict)


if __name__ == '__main__':
   a = RocAnalysis(None)
   mainid = "5eead08a17b2bf7d0fd2a202"
   # output_dir = '/mnt/ilustre/users/sanger-dev/workspace/20200803/RocAnalysis_tsg_3421_0803100328973991_2672/output/roc_analysis'
   output_dir = '/mnt/ilustre/users/sanger-dev/workspace/20200803/RocAnalysis_tsg_3421_0803100708083709_8220/output/roc_analysis'
   a.add_sg_roc(mainid, output_dir, True)
