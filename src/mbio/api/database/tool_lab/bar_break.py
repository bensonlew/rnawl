# -*- coding: utf-8 -*-
# __author__ = 'binbinzhao'

import json
from api_base import ApiBase


class BarBreak(ApiBase):
    def __init__(self, bind_object):
        super(BarBreak, self).__init__(bind_object)

    def add_bar_break(self, main_id, file_path, low_point, high_point):
        data_list_bar = []
        data_list_ishape = []
        with open(file_path) as f:
            lines = f.readlines()
            for line in lines[1:]:
                item = line.strip().split("\t")
                insert_data = {
                  "bar_break_id":self.check_objectid(main_id),
                  "name": item[0],
                  "value": float(item[1]),
                  "category": "",
                  "type": "column"
                }
                data_list_bar.append(insert_data)
                if item[2] != "NA":
                    insert_data2 = {  #  加id
                        "bar_break_id": self.check_objectid(main_id),
                        "name": item[0],
                        "group": "",
                        "type": "ishape",
                        "mean": float(item[1]),
                        "sd_high": float(item[2]),
                        "sd_low": float(item[2]),
                    }
                    data_list_ishape.append(insert_data2)
            if len(data_list_bar) == 0:
                self.bind_object.logger.info("{}文件为空！".format(file_path))
            else:
                self.col_insert_data("sg_barbreak_bar", data_list_bar)
            if len(data_list_ishape) != 0:
                self.col_insert_data("sg_barbreak_ishape", data_list_ishape)

        update_dict = {
            "column_data": json.dumps({
                "name": "name",
                "data": "value",
                "category": "category",
                "condition": {'type': "column"}
            }),
            "ishape_data": json.dumps({
                "name": "name",
                "data": ["mean", "sd_high", "sd_low"],
                "group": "group",
                "condition": {'type': "ishape"}
            }),
            "name": "column",
            "low_point": low_point,
            "high_point": high_point,
        }

        self.update_db_record("sg_barbreak", {"_id": self.check_objectid(main_id)}, update_dict)

if __name__ == '__main__':
   a = BarBreak(None)
   mainid = "5eead08a17b2bf7d0fd2d207"
   output_dir = '/mnt/ilustre/users/sanger-dev/workspace/20200618/BarBreak_tsg_3421_0618102514389705_4225/BarBreak/output/newFile.txt'
   a.add_bar_break(mainid, output_dir, 20, 30)
