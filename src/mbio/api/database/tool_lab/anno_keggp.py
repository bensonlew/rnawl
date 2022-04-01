# -*- coding: utf-8 -*-
from biocluster.config import Config
from bson import ObjectId
import json
import types
import datetime
from bson.son import SON
from api_base import ApiBase


class AnnoKeggp(ApiBase):
    def __init__(self, bind_object):
        super(AnnoKeggp, self).__init__(bind_object)


    def add_stat_detail(self, main_id, stat_path, detail_path):
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

                data = [
                    ("kegg_id", main_id),
                    ("first_category", line[0]),
                    ("second_category", line[1]),
                    ("compound_ids", line[2]),
                    ("count", int(line[3])),
                    ("type", "column")
                ]
                data = SON(data)
                data_list.append(data)
        stat_coll = self.db['metab_keggp_stat']
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

        ####detail
        #"metab_name", "compound_id", "pathway_id","description", "first_category", "second_category"
        data_list2 = list()
        with open(detail_path, "r") as f1:
            file = f1.readlines()
            for line in file[1:]:
                line = line.strip().split("\t")

                data = [
                    ("kegg_id", main_id),
                    ("name", line[0]),
                    ("compound_id", line[1]),
                    ("pathway_id", line[2]),
                    ("description", line[3]),
                    ("first_category", line[4]),
                    ("second_category", line[5])
                ]
                data = SON(data)
                data_list2.append(data)
        detail_coll = self.db['metab_keggp_detail']
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

        return('success')

    def update_main(self, main_id,):
        if not isinstance(main_id, ObjectId):
            if isinstance(main_id, types.StringTypes):
                main_id = ObjectId(main_id)
            else:
                self.bind_object.set_error('main_id必须为ObjectId对象或其对应的字符串！')
        table_tmp = [
            {"filter": "false", "field": "name", "type": "string", "sort": "false", "title": "Metabolite"},
            {"filter": "false", "field": "compound_id", "type": "string", "sort": "false", "title": "KEGG Compound ID"},
            {"filter": "false", "field": "pathway_id", "type": "string", "sort": "false", "title": "KEGG Pathway ID"},
            {"filter": "false", "field": "description", "type": "string", "sort": "false", "title": "KEGG Pathway Description"},
            {"filter": "false", "field": "first_category", "type": "string", "sort": "false", "title": "KEGG Pathway First Category"},
            {"filter": "false", "field": "second_category", "type": "string", "sort": "false", "title":"KEGG Pathway Second Category"}
        ]

        table_tmp2 = [
            {"filter": "false", "field": "first_category", "type": "string", "sort": "false", "title": "First Category"},
            {"filter": "false", "field": "second_category", "type": "string", "sort": "false", "title":"Second Category"},
            {"filter": "false", "field": "count", "type": "string", "sort": "false", "title": "Number"},
            {"filter": "false", "field": "compound_ids", "type": "string", "sort": "false", "title": ""}
        ]


        update = {
            'stat_table':json.dumps({"column": table_tmp2,"condition": {}}),
            'detail_table':json.dumps({"column": table_tmp,"condition": {}}),
            'column_data' : json.dumps({"name":"second_category","data":"count","category":"first_category", "condition":{"type":"column"}}),
            'main_id' : main_id,
        }
        self.db['metab_keggp'].update({"_id": main_id}, {"$set":update})


