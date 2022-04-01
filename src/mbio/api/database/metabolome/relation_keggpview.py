# -*- coding: utf-8 -*-
from biocluster.api.database.base import Base
from bson import ObjectId
import os
import json
import types
import datetime
import pandas as pd
from bson.son import SON

class RelationKeggpview(Base):
    def __init__(self, bind_object):
        super(RelationKeggpview, self).__init__(bind_object)
        self._project_type = "metabolome"

    def add_level(self, main_id, table):
        if not isinstance(main_id, ObjectId):
            if isinstance(main_id, types.StringTypes):
                main_id = ObjectId(main_id)
            else:
                self.bind_object.set_error('main_id必须为ObjectId对象或其对应的字符串！', code="54701602")
        data_list = list()
        self.bind_object.logger.info('开始导文件：' + table)
        with open(table, "r") as f1:
            hear = f1.readline()
            for line in f1:
                line = line.strip("\n").split("\t")
                data = [
                    ("kegg_id", main_id),
                    ("pathway_id", line[0]),
                    ("description", line[1]),
                    ("second_category", line[2]),
                    ("first_category", line[3]),
                    ("metab_list", line[4].split(';')),
                    ("metab_name", line[5].split(';')),
                    ("metab_count", int(line[6])),
                    ("gene_list", line[7].split(';')),
                    ("gene_name", line[8].split(';')),
                    ("gene_count", int(line[9])),
                    ("hyperlink", line[10])
                ]
                if len(line) == 11:
                    data.append(("pdf", line[11]))
                data = SON(data)
                data_list.append(data)
        detail_coll = self.db["relation_keggpview_level"]
        try:
            if data_list:
                detail_coll.insert_many(data_list)
            else:
                self.bind_object.logger.info("level表为空")
        except Exception,e:
            self.bind_object.set_error("导入%s信息出错： %s" , variables=(table,e), code="54701603")
        else:
            self.bind_object.logger.info("导入level信息成功")

        pathway_img = os.path.dirname(table) + '/pathway_img'
        if os.path.exists(pathway_img):
            self.add_pic(main_id, pathway_img)


    def add_pic(self, anno_kegg_id, pic_dir):
        files = os.listdir(pic_dir)
        data_list = []
        for file in files:
            if file.endswith('.html.mark'):
                pic_path = os.path.join(pic_dir, file)
            else:
                continue
            with open(pic_path, "r") as f:
                header = f.readline()
                for line in f:
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
            collection = self.db["relation_keggpview_pic"]
            if data_list:
                collection.insert_many(data_list)
        except Exception,e:
            self.bind_object.set_error("relaiton_keggpview_pic error: %s" , variables=(e))
        else:
            self.bind_object.logger.info("导relaiton_keggpview_pic成功")
