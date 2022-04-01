# -*- coding: utf-8 -*-
# __author__ = 'binbinzhao'

import json
from api_base import ApiBase


class GoBar(ApiBase):
    def __init__(self, bind_object):
        super(GoBar, self).__init__(bind_object)

    def add_sg_gobar(self, main_id, file_path):
        data_list_bar = []
        data_list_table = []
        with open(file_path) as f:
            lines = f.readlines()
            for line in lines[1:]:
                item = line.strip().split("\t")
                if ',' in item[1]:
                    name = item[1].split(',')[0]
                else:
                    name = item[1]
                insert_data = {
                    "name": name,
                    "value": item[3],
                    "category": item[2],
                    "type": "column",
                    "gobar_id": self.check_objectid(main_id),
                }
                insert_data2 = {
                    "gobar_id": self.check_objectid(main_id),
                    'go_id': item[0],
                    "description": item[1],
                    "term_type": item[2],
                    "num": item[3],
                    "percent": item[4],
                }
                data_list_bar.append(insert_data)
                data_list_table.append(insert_data2)
        if len(data_list_bar) == 0:
            self.bind_object.logger.info("{}文件为空！".format(file_path))
        else:
            self.col_insert_data("sg_gobar_bar", data_list_bar)

        if len(data_list_table) == 0:
            self.bind_object.logger.info("{}文件为空！".format(file_path))
        else:
            self.col_insert_data("sg_gobar_table", data_list_table)

        update_dict = {
            "column_data": json.dumps({
                "name": "name",
                "data": "value",
                "condition": {'type': "column"},
                "category": "category",
            }),
            "table_data": json.dumps({
                # "condition": {'type': "table"},
                "column": [
                    {"field": "go_id", "filter": "false", "sort": "false", "title": "GO ID", "type": "string"},
                    {"field": "description", "filter": "false", "sort": "false", "title": "Description",
                     "type": "string"},
                    {"field": "term_type", "filter": "true", "sort": "false", "title": "Term Type", "type": "string"},
                    {"field": "num", "filter": "false", "sort": "false", "title": "number", "type": "int"},
                    {"field": "percent", "filter": "false", "sort": "false", "title": "percent", "type": "string"}
                ],
            }),
        }
        # self.bind_object.logger.info("&&&&&&&&&&&&&&&&&&导表结束&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&")
        self.update_db_record("sg_gobar", {"_id": self.check_objectid(main_id)}, update_dict)

if __name__ == '__main__':
   a = GoBar(None)
   mainid = "5eead08a17b2bf7d0fd2d207"
   output_dir = '/mnt/ilustre/users/sanger-dev/workspace/20200706/Single_go_bar20200706093755/GoBar/output/go_table'
   a.add_sg_gobar(mainid, output_dir)
