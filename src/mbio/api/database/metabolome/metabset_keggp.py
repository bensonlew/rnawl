# -*- coding: utf-8 -*-
from biocluster.api.database.base import Base
from bson import ObjectId
import os
import json
import types
import datetime
import pandas as pd
from bson.son import SON

class MetabsetKeggp(Base):
    def __init__(self, bind_object):
        super(MetabsetKeggp, self).__init__(bind_object)
        self._project_type = "metabolome"

    def add_metabsetp(self, params, img_dir, set_name, name="KEGGP_Origin"):
        # params{set_id}用来保存代谢集主表记录
        # set_name存储代谢集名称
        # project_sn, task_id, name, status, desc, created_ts,
        task_id = self.bind_object.sheet.id
        project_sn = self.bind_object.sheet.project_sn
        created_ts = datetime.datetime.now().strftime('%Y-%m-%d %H:%M:%S')
        insert_data = {
            "project_sn": project_sn,
            "task_id": task_id,
            "desc": "代谢集KEGG通路分类",
            "created_ts": created_ts,
            "status": "end",
            "name": name,
            "params": json.dumps(params, sort_keys=True, separators=(',', ':')),
            "result_dir": img_dir,
            "set_list": set_name.split(","),
            "version" : "3.2",  #201907
            "metabset_list" : set_name  #201907
        }
        collection = self.db['metabset_keggp']
        main_id = collection.insert_one(insert_data).inserted_id
        collection.update_one({'_id': main_id}, {'$set': {'main_id': main_id}})
        return main_id

    def update_metabsetp(self, main_id, set_name):
        if not isinstance(main_id, ObjectId):
            if isinstance(main_id, types.StringTypes):
                main_id = ObjectId(main_id)
            else:
                self.bind_object.set_error('main_id必须为ObjectId对象或其对应的字符串！', code="54701601")
        collection = self.db['metabset_keggp']
        collection.update_one({'_id': main_id}, {'$set': {'set_list': set_name}})

    def add_metabsetp_level(self, main_id, table,bset_name=None):
        if not isinstance(main_id, ObjectId):
            if isinstance(main_id, types.StringTypes):
                main_id = ObjectId(main_id)
            else:
                self.bind_object.set_error('main_id必须为ObjectId对象或其对应的字符串！', code="54701602")
        data_list = list()
        self.bind_object.logger.info('开始导文件：' + table)
        with open(table, "r") as f1:
            file = f1.readlines()
            for line in file[1:]:
                line = line.strip().split("\t")
                data = [
                    ("kegg_id", main_id),
                    ("pathway_id", line[0]),
                    ("description", line[1]),
                    ("first_category", line[2]),
                    ("second_category", line[3]),
                    ("compound_id", line[4]),
                    ("metab_id", line[5]),
                    ("hyperlink", line[6]),
                    ("count", int(float(line[7])))
                ]
                if len(line)>9:
                    if bset_name:
                        if len(line) < 11:
                            self.bind_object.logger.info('缺少列')
                            self.bind_object.logger.info(line)
                        data.extend([
                            ("compound_id_"+bset_name, line[8]),
                            ("metab_id_"+bset_name, line[9]),
                            ("count_"+bset_name, int(float(line[10])))
                    ])

                data = SON(data)
                data_list.append(data)
        detail_coll = self.db["metabset_keggp_level"]
        try:
            if data_list:
                detail_coll.insert_many(data_list)
            else:
                #self.db['metabset_keggp'].remove({"_id": main_id})
                self.bind_object.logger.info("level表为空")
        except Exception,e:
            self.bind_object.set_error("导入%s信息出错： %s" , variables=(table,e), code="54701603")
        else:
            self.bind_object.logger.info("导入level信息成功")

    def add_metabsetp_stat(self, main_id, set_name, table):
        if not isinstance(main_id, ObjectId):
            if isinstance(main_id, types.StringTypes):
                main_id = ObjectId(main_id)
            else:
                self.bind_object.set_error('main_id必须为ObjectId对象或其对应的字符串！', code="54701604")
        data_list = list()
        with open(table, "r") as f1:
            self.bind_object.logger.info("start import star table1")
            file = f1.readlines()
            self.bind_object.logger.info("start import star table2")
            for line in file[1:]:
                self.bind_object.logger.info("start import star table3")
                line = line.strip().split("\t")
                data = [
                    ("kegg_id", main_id),
                    ("first_category", line[0]),
                    ("second_category", line[1]),
                    ("metab_id", line[2]),
                    ("count", int(float(line[3]))),
                    ("set_name", set_name)
                ]
                data = SON(data)
                data_list.append(data)
        detail_coll = self.db["metabset_keggp_stat"]
        self.bind_object.logger.info("start import star table")
        try:
            if data_list:
                detail_coll.insert_many(data_list)
            else:
                #self.db['metabset_keggp'].remove({"_id": main_id})
                self.bind_object.logger.info("stat 表为空")
        except Exception,e:
            self.bind_object.set_error("导入%s信息出错： %s" , variables=(table,e), code="54701605")
        else:
            self.bind_object.logger.info("导入stat信息成功")

        ##导 metabset_keggp_pic
        pathway_img = os.path.dirname(table) + '/pathway_img'
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
                    if len(line) < 6:
                        self.bind_object.logger.info(line)
                    data = [
                        ("kegg_id", anno_kegg_id),
                        ("pathway_id", line[0]),
                        ("query", line[1]),
                        ("shape", line[2]),
                        ("coords", line[3]),
                        ("title", line[4]),
                        ("href", line[5]),

                    ]
                    if len(line) >6:
                        data.append(("color", line[6]))
                    data = SON(data)
                    data_list.append(data)
        try:
            collection = self.db["metabset_keggp_pic"]
            if data_list:
                collection.insert_many(data_list)
        except Exception,e:
            self.bind_object.set_error("metabset_keggp_pic error: %s" , variables=(e))
        else:
            self.bind_object.logger.info("导metabset_keggp_pic成功")