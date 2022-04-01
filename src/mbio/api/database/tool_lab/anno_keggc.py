# -*- coding: utf-8 -*-
from biocluster.config import Config
from bson import ObjectId
import json
import types
import datetime
from bson.son import SON
from api_base import ApiBase


class AnnoKeggc(ApiBase):
    def __init__(self, bind_object):
        super(AnnoKeggc, self).__init__(bind_object)


    def add_stat_detail(self, main_id, stat_path, sub_database):
        if not isinstance(main_id, ObjectId):
            if isinstance(main_id, types.StringTypes):
                main_id = ObjectId(main_id)
            else:
                self.bind_object.set_error('main_id必须为ObjectId对象或其对应的字符串！')
        data_list = list()
        with open(stat_path, "r") as f1:
            file = f1.readlines()
            for line in file[1:]:
                line = line.strip().split("\t")
                second_c = line[1]
                if second_c == 'Others':
                    second_c = line[0] + '_Others'
                data = [
                    ("kegg_id", main_id),
                    ("first_category", line[0]),
                    ("second_category", second_c),
                    ("hyperlink", line[2]),
                    ("names", line[3]),
                    ("count", int(line[4])),
                    ("database", sub_database),
                    ("type", "column")
                ]
                data = SON(data)
                data_list.append(data)
        stat_coll = self.db['keggc_detail']
        try:
            if data_list:
                stat_coll.insert_many(data_list)
            else:
                self.bind_object.logger.info("stat表为空")
                return('no data')
        except Exception as e:
            self.bind_object.set_error("导入%s信息出错： %s" , variables=(stat_path,e))
        else:
            self.bind_object.logger.info("导入stat信息成功")
        return('success')

    def update_main(self, main_id, true_database):
        if not isinstance(main_id, ObjectId):
            if isinstance(main_id, types.StringTypes):
                main_id = ObjectId(main_id)
            else:
                self.bind_object.set_error('main_id必须为ObjectId对象或其对应的字符串！')
        table_tmp = [
                    {"filter": "false", "field": "first_category", "type": "string", "sort": "false", "title": "First Category"},
                    {"filter": "false", "field": "second_category", "type": "string", "sort": "false", "title":"Second Category"},
                    {"filter": "false", "field": "count", "type": "string", "sort": "false", "title": "Number"}
                ]
        update = {
            'table_data':json.dumps({"column": table_tmp,"condition": {}}),
            'column_data' : json.dumps({"name":"second_category","data":"count","category":"first_category", "condition":{"type": "column"}}),
            'column_names' : json.dumps({"data":true_database}),
            'main_id' : main_id,
        }
        self.db['keggc'].update({"_id": main_id}, {"$set":update})


