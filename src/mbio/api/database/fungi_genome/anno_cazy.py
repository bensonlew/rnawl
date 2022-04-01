# -*- coding: utf-8 -*-
# __author__ = 'zhujuan'
# last modify: 2018.03.10

from biocluster.api.database.base import Base, report_check
import datetime
import json
from bson.son import SON
import pandas as pd
from bson.objectid import ObjectId

class AnnoCazy(Base):
    def __init__(self, bind_object):
        super(AnnoCazy, self).__init__(bind_object)
        self._project_type = "fungigenome"

    @report_check
    def add_anno_cazy(self, params=None):
        task_id = self.bind_object.sheet.id
        project_sn = self.bind_object.sheet.project_sn
        created_ts = datetime.datetime.now().strftime('%Y-%m-%d %H:%M:%S')
        insert_data = {
            "project_sn": project_sn,
            "task_id": task_id,
            "status": "end",
            "desc": "CAZy注释",
            "params": json.dumps(params, sort_keys=True, separators=(',', ':')),
            "created_ts": created_ts,
            "name": "AnnoCAZy_Origin",
            "settled_params": json.dumps({"version": "cazy_v20200408"})
        }
        collection = self.db["anno_cazy"]
        inserted_id = collection.insert_one(insert_data).inserted_id
        return inserted_id

    @report_check
    def add_anno_cazy_detail(self, inserted_id, specimen_id, anno):
        data_list = []
        ann = pd.read_table(anno, sep='\t', header=0)
        if len(ann) < 1:
            return
        for i in range(len(ann)):
            data = {
                "cazy_id": ObjectId(inserted_id),
                "specimen_id": specimen_id,
                "gene_id":  ann["Gene ID"][i],
                "location": ann["Location"][i],
                "family": ann["Family"][i],
                "class": ann["Class"][i],
                "evalue": ann["Evalue"][i],
            }
            if "Coverd_fraction" in ann.columns:
                data["covered"] = ann["Coverd_fraction"][i]
            elif "Coverage(%)" in ann.columns :
                data["covered"] = ann["Coverage(%)"][i]
            data_son = SON(data)
            data_list.append(data_son)
        try:
            collection = self.db["anno_cazy_detail"]
            collection.insert_many(data_list)

            main_collection = self.db['anno_cazy']
            main_collection.update({"_id": ObjectId(inserted_id)}, {"$set": {"status": 'end',"origin_id": ObjectId(inserted_id)}})
        except Exception, e:
            self.bind_object.logger.info("导入%s结果表出错:%s" % (anno, e))
        else:
            self.bind_object.logger.info("导入%s结果表成功!" % anno)
