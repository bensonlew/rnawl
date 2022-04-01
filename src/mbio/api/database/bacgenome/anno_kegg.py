# -*- coding: utf-8 -*-
# __author__ = 'zhujuan'
# last modify: 2018.03.10

from biocluster.api.database.base import Base, report_check
import datetime
import json
from bson.son import SON
import pandas as pd
from bson.objectid import ObjectId
import os
import re

class AnnoKegg(Base):
    def __init__(self, bind_object):
        super(AnnoKegg, self).__init__(bind_object)
        self._project_type = "bacgenome"

    @report_check
    def add_anno_kegg(self, path1, path2, path3, params=None):
        task_id = self.bind_object.sheet.id
        project_sn = self.bind_object.sheet.project_sn
        created_ts = datetime.datetime.now().strftime('%Y-%m-%d %H:%M:%S')
        insert_data = {
            "project_sn": project_sn,
            "task_id": task_id,
            "status": "end",
            "desc": "KEGG注释",
            "params": json.dumps(params, sort_keys=True, separators=(',', ':')),
            "created_ts": created_ts,
            "file_path": [path1, path2, path3],
            "name": "AnnoKEGG_Origin",
            "version" : "3.1",
            "settled_params": json.dumps({"version": "kegg_v94.2"})
        }
        collection = self.db["anno_kegg"]
        inserted_id = collection.insert_one(insert_data).inserted_id
        collection.update({'_id': inserted_id},{'$set':{'main_id':inserted_id}}) #guanqing 20180813
        return inserted_id

    @report_check
    def add_anno_kegg_detail(self, inserted_id, specimen_id, anno):
        data_list = []
        ann = pd.read_table(anno, sep='\t', header=0)
        if len(ann) < 1:
            return
        for i in range(len(ann)):
            data = {
                "kegg_id": ObjectId(inserted_id),
                "specimen_id": specimen_id,
                "hit_gene":ann["Kegg_name_id"][i],
                "gene_id": ann["Gene ID"][i],
                "location": ann["Location"][i],
                "ko": ann["KO"][i],
                "gene_name": ann["Gene"][i],
                "ko_des": ann["Definition"][i],
                "pathway": ann["Pathway"][i],
                "enzyme": ann["Enzyme"][i],
                "en_des": ann["Enzyme_description"][i],
                "pathway_des": ann["Level3"][i],
            }
            data_son = SON(data)
            data_list.append(data_son)
        try:
            collection = self.db["anno_kegg_detail"]
            collection.insert_many(data_list)
        except Exception, e:
            self.bind_object.logger.info("导入%s结果表出错:%s" % (anno, e))
        else:
            self.bind_object.logger.info("导入%s结果表成功!" % anno)

    @report_check
    def add_anno_kegg_level(self, inserted_id, specimen_id, level):
        des_pat = re.compile('(C\d*)\s*\((.*)\)')
        data_list = []
        compound_table = self.db['anno_kegg_pic']
        ann = pd.read_table(level, sep='\t', header=0)
        if len(ann) < 1:
            return
        for i in range(len(ann)):
            if ann["Level1"][i] != "-":
                data = {
                    "kegg_id": ObjectId(inserted_id),
                    "specimen_id": specimen_id,
                    "level1": ann["Level1"][i],
                    "level2": ann["Level2"][i],
                    "level3": ann["Level3"][i],
                    "gene_num": ann["Gene nu"][i],
                    "gene_list": ann["Gene list"][i],
                    "pathway": ann["Pathway list"][i],
                    "ko": ann["KO list"][i]
                }
                pathway = ann["Pathway list"][i].split(';')[0]
                res = compound_table.find({"kegg_id": ObjectId(inserted_id),"specimen_id": specimen_id,"pathway_id":pathway,"shape" : "circle"})
                c_des_list = []
                c_id_list = []
                if res:
                    for rs in res:
                        f = des_pat.findall(rs["title"])
                        if f:
                            c_id = f[0][0]
                            c_des = f[0][1]
                            c_des_list.append(c_des)
                            c_id_list.append(c_id)
                data['compound'] = ';'.join(c_id_list)
                data['comp_des'] = ';'.join(c_des_list)

                data_son = SON(data)
                data_list.append(data_son)
        try:
            collection = self.db["anno_kegg_level"]
            collection.insert_many(data_list)
        except Exception, e:
            self.bind_object.logger.info("导入%s结果表出错:%s" % (level, e))
        else:
            self.bind_object.logger.info("导入%s结果表成功!" % level)

    @report_check
    def add_anno_kegg_pic(self, anno_kegg_id,specimen_id, pic_dir):    #guanqing.zou 20190404
        files = os.listdir(pic_dir)
        data_list = []
        for file in files:
            if file.endswith('.html.mark'):
                pic_path = os.path.join(pic_dir, file)
            else:
                continue
            with open(pic_path, "r") as f:
                lines = f.readlines()
                for line in lines:
                    line = line.strip().split("\t")
                    if line[1]=='poly' and line[2]=='':  # 没有注释到的K的直线信息不导数据库。减少数据库无用信息的存储
                        continue

                    data = [
                        ("kegg_id", anno_kegg_id),
                        ("pathway_id", line[0]),
                        ("shape", line[1]),
                        ("color", line[2]),
                        ("coords", line[4]),
                        ("title", line[5]),
                        ("query", line[6]),
                        ("specimen_id",specimen_id)
                        #("href", line[7])
                    ]
                    data = SON(data)
                    data_list.append(data)
        try:
            collection = self.db["anno_kegg_pic"]
            collection.insert_many(data_list)
        except Exception,e:
            self.bind_object.set_error("anno_kegg_pic error: %s" , variables=(e))
        else:
            self.bind_object.logger.info("导入anno_kegg_pic成功")