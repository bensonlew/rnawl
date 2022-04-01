# -*- coding: utf-8 -*-
from biocluster.api.database.base import Base
from bson import ObjectId
import pandas as pd
import os
import json
import types
import datetime
from bson.son import SON


class AnnoKeggp(Base):
    def __init__(self, bind_object):
        super(AnnoKeggp, self).__init__(bind_object)
        self._project_type = 'metabolome'

    def export_test(self, table_id, table_type, organism, level_path, stat_path):
        params = {
            "table_id": str(table_id),
            "table_type": table_type,
            "organism": organism
        }
        kegg_main = self.add_anno_keggp_main("KEGGP_", "/mnt/ilustre/users/sanger-dev/sg-users/guhaidong/Metabolite/Annop/pathway_img", params)
        self.add_level_detail(kegg_main, level_path)
        self.add_stat_detail(kegg_main, stat_path)

    def add_anno_keggp_main(self, name, result_dir, params):
        import datetime
        created_ts = datetime.datetime.now().strftime('%Y-%m-%d %H:%M:%S')
        insert_data = {
            "project_sn": self.bind_object.sheet.project_sn,
            "task_id": self.bind_object.sheet.id,
            "desc": "kegg pathway注释表",
            "name": name,
            "created_ts": created_ts,
            "params": json.dumps(params, sort_keys=True, separators=(',', ':')),
            "status": "end",
            "result_dir": result_dir,
            "version" : "3.2"
        }
        collection = self.db['anno_keggp']
        main_id = collection.insert_one(insert_data).inserted_id
        collection.update_one({'_id': main_id}, {'$set': {'main_id': main_id}})
        return main_id

    def add_level_detail(self, main_id, level_path):
        if not isinstance(main_id, ObjectId):
            if isinstance(main_id, types.StringTypes):
                main_id = ObjectId(main_id)
            else:
                self.bind_object.set_error("main_id必须为ObjectId对象或其他对应的字符串！", code="54701801")
        data_list = list()
        tmp = level_path.split("/")
        img_dir = "/".join(tmp[0:len(tmp)-1])
        img_path = os.path.join(img_dir, "pathway_img")
        with open(level_path, "r") as f1:
            file = f1.readlines()
            for line in file[1:]:
                line = line.strip().split("\t")
                if line[2] == "Global and overview maps":
                    type = "global"
                else:
                    type = "important"
                pdf = os.path.join(img_path,line[0].replace("map","ko")+".pdf")
                self.bind_object.logger.info(pdf)
                if os.path.exists(pdf):
                    pdf = line[0]
                else:
                    pdf = "-"
                data = [
                    ('kegg_id', main_id),
                    ('pathway_id', line[0]),
                    ('description', line[3]),
                    ('first_category', line[1]),
                    ('second_category', line[2]),
                    ('compound_id', line[4]),
                    ('metab_id', line[5]),
                    ('count', int(line[6])),
                    ('hyperlink', line[7]),
                    ('type', type),
                    ("pdf", pdf)
                ]
                data = SON(data)
                data_list.append(data)

        level_collection = self.db['anno_keggp_level']
        try:
            if data_list:
                level_collection.insert_many(data_list)
        except Exception, e:
            self.bind_object.set_error("keggp level error: %s" , variables=(e), code="54701802")
        else:
            self.bind_object.logger.info("anno_keggp_level导表成功")

    def add_stat_detail(self, main_id, stat_path):
        if not isinstance(main_id, ObjectId):
            if isinstance(main_id, types.StringTypes):
                main_id = ObjectId(main_id)
            else:
                self.bind_object.set_error('main_id必须为ObjectId对象或其对应的字符串！', code="54701803")
        data_list = list()
        with open(stat_path, "r") as f1:
            file = f1.readlines()
            for line in file[1:]:
                line = line.strip().split("\t")
                data = [
                    ("kegg_id", main_id),
                    ("first_category", line[0]),
                    ("second_category", line[1]),
                    ("metab_id", line[2]),
                    ("count", line[3]),

                ]
                data = SON(data)
                data_list.append(data)
        stat_collection = self.db['anno_keggp_stat']
        try:
            if data_list:
                stat_collection.insert_many(data_list)
        except Exception, e:
            self.bind_object.set_error("keggp stat error: %s" , variables=(e), code="54701804")
        else:
            self.bind_object.logger.info("anno_keggp_stat 成功")
        ## 导kegg_pic
        pathway_img = os.path.dirname(stat_path) + '/pathway_img'
        if os.path.exists(pathway_img):
            self.add_anno_kegg_pic(main_id, pathway_img)


    def add_anno_kegg_pic(self, anno_kegg_id, pic_dir):    #guanqing.zou 20190531
        files = os.listdir(pic_dir)
        data_list = []
        for file in files:
            if file.endswith('.html.mark'):
                pic_path = os.path.join(pic_dir, file)
            else:
                continue
            with open(pic_path, "r") as f:
                lines = f.readlines()
                for line in lines[1:]:
                    line = line.strip().split("\t")
                    data = [
                        ("kegg_id", anno_kegg_id),
                        ("pathway_id", line[0]),
                        ("query", line[1]),
                        ("shape", line[2]),
                        #("color", line[2]),
                        ("coords", line[3]),
                        ("title", line[4]),
                        ("href", line[5])
                    ]
                    data = SON(data)
                    data_list.append(data)
        try:
            collection = self.db["anno_keggp_pic"]
            if data_list:
                collection.insert_many(data_list)
        except Exception,e:
            self.bind_object.set_error("anno_keggp_pic error: %s" , variables=(e))
        else:
            self.bind_object.logger.info("导入anno_keggp_pic成功")


