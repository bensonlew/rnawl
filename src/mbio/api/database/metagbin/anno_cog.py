# -*- coding: utf-8 -*-
# __author__ = 'gaohao'

from biocluster.api.database.base import Base, report_check
import datetime
import json
from bson.son import SON
import pandas as pd
from bson.objectid import ObjectId


class AnnoCog(Base):
    def __init__(self, bind_object):
        super(AnnoCog, self).__init__(bind_object)
        self._project_type = "metagbin"

    @report_check
    def add_anno_cog(self, params=None, name=None, task_id=None):
        if task_id is None:
            task_id = self.bind_object.sheet.id
        project_sn = self.bind_object.sheet.project_sn
        created_ts = datetime.datetime.now().strftime('%Y-%m-%d %H:%M:%S')
        insert_data = {
            "project_sn": project_sn,
            "task_id": task_id,
            "status": "end",
            "desc": "COG注释",
            "params": json.dumps(params, sort_keys=True, separators=(',', ':')),
            "created_ts": created_ts,
            "name": name if name else "AnnoCOG_Origin"
        }
        collection = self.db["anno_cog"]
        inserted_id = collection.insert_one(insert_data).inserted_id
        collection.update({'_id': inserted_id},{'$set':{'main_id': inserted_id}})
        return inserted_id

    @report_check
    def add_anno_cog_detail(self, inserted_id, genome_id, anno):
        data_list = []
        ann = pd.read_table(anno, sep='\t', header=0)
        if len(ann) < 1:
            return
        for i in range(len(ann)):
            data = {
                "cog_id": ObjectId(inserted_id),
                "genome_id": genome_id,
                "gene_id":  ann["Gene ID"][i],
                "location": ann["Location"][i],
                "gene_des": ann["Gene Description"][i],
                "cog": ann["COG ID"][i],
                "cog_des": ann["COG Description"][i],
                "type": ann["COG Type"][i],
                "type_dec": ann["Type Description"][i],
                "category": ann["Category"][i],
            }
            data_son = SON(data)
            data_list.append(data_son)
        try:
            collection = self.db["anno_cog_detail"]
            collection.insert_many(data_list)
        except Exception, e:
            self.bind_object.logger.info("导入%s结果表出错:%s" % (anno, e))
        else:
            self.bind_object.logger.info("导入%s结果表成功!" % anno)
