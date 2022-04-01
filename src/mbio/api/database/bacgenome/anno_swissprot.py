# -*- coding: utf-8 -*-
# __author__ = 'zhujuan'
# last modify: 2018.03.10

from biocluster.api.database.base import Base, report_check
import datetime
import json
from bson.son import SON
import pandas as pd
from bson.objectid import ObjectId

class AnnoSwissprot(Base):
    def __init__(self, bind_object):
        super(AnnoSwissprot, self).__init__(bind_object)
        self._project_type = "bacgenome"

    @report_check
    def add_anno_swissprot(self, params=None):
        task_id = self.bind_object.sheet.id
        project_sn = self.bind_object.sheet.project_sn
        created_ts = datetime.datetime.now().strftime('%Y-%m-%d %H:%M:%S')
        insert_data = {
            "project_sn": project_sn,
            "task_id": task_id,
            "status": "end",
            "desc": "Swiss-prot注释",
            "params": json.dumps(params, sort_keys=True, separators=(',', ':')),
            "created_ts": created_ts,
            "name": "AnnoSwissprot_Origin",
            "settled_params": json.dumps({"version": "swissprot_v20200617"})
        }
        collection = self.db["anno_swissprot"]
        inserted_id = collection.insert_one(insert_data).inserted_id
        collection.update({'_id': inserted_id},{'$set':{'main_id':inserted_id}}) #guanqing 20180813
        return inserted_id

    @report_check
    def add_anno_swissprot_detail(self, inserted_id, specimen_id, anno):
        data_list = []
        ann = pd.read_table(anno, sep='\t', header=0)
        if len(ann) < 1:
            return
        for i in range(len(ann)):
            data = {
                "swissprot_id": ObjectId(inserted_id),
                "specimen_id": specimen_id,
                "gene_id":  ann["Gene ID"][i],
                "location": ann["Location"][i],
                "swiss_hit": ann["Hit"][i],
                "swiss_des": ann["Hit-Description"][i],
                "len": int(float(ann["Gene Len"][i])+1)*3 if ann["Gene Len"][i] != "-" else "-",
                "start": int(float(ann["Gene Start"][i])) if ann["Gene Start"][i] != "-" else "-",
                "end": int(float(ann["Gene End"][i])) if ann["Gene End"][i] != "-" else "-",
                "s_start": int(float(ann["Hit Start"][i])) if ann["Hit Start"][i] != "-" else "-",
                "s_end": int(float(ann["Hit End"][i])) if ann["Hit End"][i] != "-" else "-",
                "hit_len": int(float(ann["Hit Len"][i])) if ann["Hit Len"][i] != "-" else "-",
                "identity": float(ann["Identity"][i]) if ann["Identity"][i] != "-" else "-",
                "evalue": float(ann["Evalue"][i]) if ann["Evalue"][i] != "-" else "-",
                "score": float(ann["Score"][i]) if ann["Score"][i] != "-" else "-",
                "coverage" : round(abs(float(ann["Gene Start"][i]) - float(ann["Gene End"][i]))/float(ann["Gene Len"][i]) ,4)*100 if ann["Gene Len"][i] != "-" else "-",
            }
            data_son = SON(data)
            data_list.append(data_son)
        try:
            collection = self.db["anno_swissprot_detail"]
            collection.insert_many(data_list)
        except Exception, e:
            self.bind_object.logger.info("导入%s结果表出错:%s" % (anno, e))
        else:
            self.bind_object.logger.info("导入%s结果表成功!" % anno)
