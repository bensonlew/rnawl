# -*- coding: utf-8 -*-
# __author__ = 'zhujuan'
# last modify: 2018.03.10

from biocluster.api.database.base import Base, report_check
import datetime
import json
from bson.son import SON
import pandas as pd
from bson.objectid import ObjectId

class AnnoNr(Base):
    def __init__(self, bind_object):
        super(AnnoNr, self).__init__(bind_object)
        self._project_type = "fungigenome"

    @report_check
    def add_anno_nr(self, params=None):
        task_id = self.bind_object.sheet.id
        project_sn = self.bind_object.sheet.project_sn
        created_ts = datetime.datetime.now().strftime('%Y-%m-%d %H:%M:%S')
        insert_data = {
            "project_sn": project_sn,
            "task_id": task_id,
            "status": "end",
            "desc": "NR注释",
            "params": json.dumps(params, sort_keys=True, separators=(',', ':')),
            "created_ts": created_ts,
            "name": "AnnoNR_Origin",
            "settled_params": json.dumps({"version": "nr_v20200604"})
        }
        collection = self.db["anno_nr"]
        inserted_id = collection.insert_one(insert_data).inserted_id
        return inserted_id

    @report_check
    def add_anno_nr_detail(self, inserted_id, specimen_id, anno):
        data_list = []
        ann = pd.read_table(anno, sep='\t', header=0)
        if len(ann) < 1:
            return
        for i in range(len(ann)):
            data = {
                "nr_id": ObjectId(inserted_id),
                "specimen_id": specimen_id,
                "gene_id":  ann["Gene ID"][i],
                "location": ann["Location"][i],
                "nr_hit": ann["Hit"][i],
                "gene_des": ann["Hit-Description"][i],
                "len": int(float(ann["Gene Len"][i])) if ann["Gene Len"][i] != "-" else "-",
                "start": int(float(ann["Gene Start"][i])) if ann["Gene Start"][i] != "-" else "-",
                "end": int(float(ann["Gene End"][i])) if ann["Gene End"][i] != "-" else "-",
                "s_start": int(float(ann["Hit Start"][i])) if ann["Hit Start"][i] != "-" else "-",
                "s_end": int(float(ann["Hit End"][i])) if ann["Hit End"][i] != "-" else "-",
                "hit_len": int(float(ann["Hit Len"][i])) if ann["Hit Len"][i] != "-" else "-",
                "identity": float(ann["Identity"][i]) if ann["Identity"][i] != "-" else "-",
                "evalue": float(ann["Evalue"][i]) if ann["Evalue"][i] != "-" else "-",
                "score": float(ann["Score"][i]) if ann["Score"][i] != "-" else "-",
            }
            data_son = SON(data)
            data_list.append(data_son)
        try:
            collection = self.db["anno_nr_detail"]
            collection.insert_many(data_list)
            main_collection = self.db["anno_nr"]
            main_collection.update({"_id": ObjectId(inserted_id)},
                                   {"$set": {"status": 'end',"origin_id": ObjectId(inserted_id)}})
        except Exception, e:
            self.bind_object.logger.info("导入%s结果表出错:%s" % (anno, e))
        else:
            self.bind_object.logger.info("导入%s结果表成功!" % anno)
