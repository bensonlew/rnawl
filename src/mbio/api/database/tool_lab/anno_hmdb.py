# -*- coding: utf-8 -*-
from biocluster.config import Config
from bson import ObjectId
import json
import types
import datetime
from bson.son import SON
from api_base import ApiBase


class AnnoHmdb(ApiBase):
    def __init__(self, bind_object):
        super(AnnoHmdb, self).__init__(bind_object)

    def add_detail(self,main_id, detail_path):
        if not isinstance(main_id, ObjectId):
            if isinstance(main_id, types.StringTypes):
                main_id = ObjectId(main_id)
            else:
                self.bind_object.set_error('main_id必须为ObjectId对象或其对应的字符串！')

        ####detail
        #["metab_id","Metabolite","HMDB ID", "Kingdom", "Superclass", "Class", "Subclass"]
        data_list2 = list()
        with open(detail_path, "r") as f1:
            file = f1.readlines()
            for line in file[1:]:
                line = line.strip().split("\t")

                data = [
                    ("anno_id", main_id),
                    ("metab_id", line[0]),
                    ("name", line[1]),
                    ("hmdb", line[2]),
                    ("kingdom", line[3]),
                    ("superclass", line[4]),
                    ("class", line[5]),
                    ("subclass", line[6])
                ]
                data = SON(data)
                data_list2.append(data)
        detail_coll = self.db['metab_hmdb_detail']
        try:
            if data_list2:
                detail_coll.insert_many(data_list2)
            else:
                self.bind_object.logger.info("detail表为空")
                return('no data')
        except Exception as e:
            self.bind_object.set_error("导入%s信息出错： %s" , variables=(detail_path,e))
        else:
            self.bind_object.logger.info("导入detail信息成功")


    def add_stat(self, main_id, stat_path, level):
        if not isinstance(main_id, ObjectId):
            if isinstance(main_id, types.StringTypes):
                main_id = ObjectId(main_id)
            else:
                self.bind_object.set_error('main_id必须为ObjectId对象或其对应的字符串！')

        data_list = list()
        #[levelname, "Number", "Metab_ids"]
        with open(stat_path, "r") as f1:
            file = f1.readlines()
            for id,line in enumerate(file[1:],1):
                line = line.strip().split("\t")

                data = [
                    ("anno_id", main_id),
                    ("name", line[0]),
                    ("number", int(line[1])),
                    ("metab_ids", line[2]),
                    ("level",level),
                    ("rank",  id)
                ]
                data = SON(data)
                data_list.append(data)
        stat_coll = self.db['metab_hmdb_stat']
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

    def update_main(self, main_id,):
        if not isinstance(main_id, ObjectId):
            if isinstance(main_id, types.StringTypes):
                main_id = ObjectId(main_id)
            else:
                self.bind_object.set_error('main_id必须为ObjectId对象或其对应的字符串！')
        table_tmp = [
            {"filter": "false", "field": "name", "type": "string", "sort": "false", "title": "Metabolite"},
            {"filter": "false", "field": "hmdb", "type": "string", "sort": "false", "title": "HMDB ID"},
            {"filter": "false", "field": "kingdom", "type": "string", "sort": "false", "title": "Kingdom"},
            {"filter": "false", "field": "superclass", "type": "string", "sort": "false", "title": "HMDB SuperClass"},
            {"filter": "false", "field": "class", "type": "string", "sort": "false", "title": "HMDB Class"},
            {"filter": "false", "field": "subclass", "type": "string", "sort": "false", "title":"HMDB SubClass"}
        ]

        table_tmp2 = [
            {"filter": "false", "field": "name", "type": "string", "sort": "false", "title": "Name"},
            {"filter": "false", "field": "number", "type": "string", "sort": "false", "title": "Number"},
            {"filter": "false", "field": "metab_ids", "type": "string", "sort": "false", "title": ""},
        ]


        update = {
            'detail_table':json.dumps({"column": table_tmp,"condition": {}}),
            'stat_table':json.dumps({"column": table_tmp2,"condition": {}}),
            'pie_data' : json.dumps({"name":"name","data":"number", "rank":"rank","condition":{}}),
            'col_name' :  json.dumps({"data": ["superclass", "class", "subclass"]}),
            'main_id' : main_id,
        }
        self.db['metab_hmdb'].update({"_id": main_id}, {"$set":update})


