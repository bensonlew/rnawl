# -*- coding: utf-8 -*-
# __author__ = 'zhujuan'
# last modify: 2018.03.10

from biocluster.api.database.base import Base, report_check
import datetime
import json
from bson.son import SON
import pandas as pd
from bson.objectid import ObjectId


class AnnoKegg(Base):
    def __init__(self, bind_object):
        super(AnnoKegg, self).__init__(bind_object)
        self._project_type = "fungigenome"

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
            "settled_params": json.dumps({"version": "kegg_v94.2"})
        }
        collection = self.db["anno_kegg"]
        inserted_id = collection.insert_one(insert_data).inserted_id
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
                "gene_id": ann["Gene ID"][i],
                "location": ann["Location"][i],
                "ko": ann["KO"][i],
                "gene_name": ann["Gene"][i],
                "ko_des": ann["Definition"][i],
                "pathway": ann["Pathway"][i],
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
        data_list = []
        ann = pd.read_table(level, sep='\t', header=0)
        if len(ann) < 1:
            return
        for i in range(len(ann)):
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
            data_son = SON(data)
            data_list.append(data_son)
        try:
            collection = self.db["anno_kegg_level"]
            collection.insert_many(data_list)
            main_collection = self.db["anno_kegg"]
            main_collection.update({"_id": ObjectId(inserted_id)},
                                   {"$set": {"status": 'end',"origin_id": ObjectId(inserted_id)}})
        except Exception, e:
            self.bind_object.logger.info("导入%s结果表出错:%s" % (level, e))
        else:
            self.bind_object.logger.info("导入%s结果表成功!" % level)
