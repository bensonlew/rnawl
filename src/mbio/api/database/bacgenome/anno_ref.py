# -*- coding: utf-8 -*-
# __author__ = 'zhujuan'
# last modify: 2018.03.10

from biocluster.api.database.base import Base, report_check
import datetime
import json
from bson.son import SON
import pandas as pd
import os, re
from bson.objectid import ObjectId


class AnnoRef(Base):
    def __init__(self, bind_object):
        super(AnnoRef, self).__init__(bind_object)
        self._project_type = "bacgenome"

    @report_check
    def add_anno_ref(self, specimen_id, params=None):
        task_id = self.bind_object.sheet.id
        project_sn = self.bind_object.sheet.project_sn
        created_ts = datetime.datetime.now().strftime('%Y-%m-%d %H:%M:%S')
        insert_data = {
            "project_sn": project_sn,
            "task_id": task_id,
            "status": "end",
            "desc": "参考基因组注释",
            "params": json.dumps(params, sort_keys=True, separators=(',', ':')),
            "created_ts": created_ts,
            "specimen_id": specimen_id
        }
        collection = self.db["anno_ref"]
        inserted_id = collection.insert_one(insert_data).inserted_id
        collection.update({'_id': inserted_id},{'$set':{'main_id':inserted_id}}) #guanqing 20180813
        return inserted_id

    @report_check
    def add_anno_ref_detail(self, inserted_id, specimen_id, anno, task_id):
        data_list = []
        ann = pd.read_table(anno, sep='\t', header=0)
        name = re.subn("\.", "_", os.path.basename(anno).rsplit(".", 1)[0].rsplit("_")[-1])[0]
        ref = "ref_" + name
        refdes = "refdes_" + name
        summary_collection = self.db["anno_summary"]
        result = summary_collection.find_one({"task_id": task_id})
        summary_id = result['_id']
        if name not in result['ref']:
            summary_collection.update({"task_id": task_id}, {"$push": {"ref": name}}, False, True)
        summary_detail_collection = self.db["anno_summary_detail"]
        geneid = {}
        for i in range(len(ann)):
            data = {
                "ref_id": ObjectId(inserted_id),
                "gene_id": ann["Gene ID"][i],
                # "location": ann["Location"][i],
                "ref": ann["Hit"][i],
                "refdes": ann["Hit-Description"][i],
                "len": ann["Gene Len"][i],
                "start": ann["Gene Start"][i],
                "end": ann["Gene End"][i],
                "s_len": ann["Ref Len"][i],
                "s_start": ann["Hit Start"][i],
                "s_end": ann["Hit End"][i],
                "hit_len": ann["Hit Len"][i],
                "identity": ann["Identity"][i],
                "evalue": ann["Evalue"][i],
                "score": ann["Score"][i],
                "coverage": round(abs(float(ann["Gene End"][i])-float(ann["Gene Start"][i]))/float(ann["Gene Len"][i]),4)*100
            }
            if ann["Gene ID"][i] not in geneid.keys():
                summary_detail_collection.update_one(
                    {"summary_id": ObjectId(summary_id), "specimen_id": specimen_id, "gene_id": ann["Gene ID"][i]},
                    {"$set": {ref: ann["Hit"][i], refdes: ann["Hit-Description"][i]}})
                geneid[ann["Gene ID"][i]] = 1
            data_son = SON(data)
            data_list.append(data_son)
        try:
            collection = self.db["anno_ref_detail"]
            collection.insert_many(data_list)
            main_col = self.db["anno_ref"]
            main_col.update_one({"_id":ObjectId(inserted_id)},{"$set":{"main_id":ObjectId(inserted_id)}})  #update main_id
            summary_stat = self.db["anno_summary_stat"]
            hit_ref_nums = ann["Gene ID"].drop_duplicates().count()
            summary_stat.update_one({"summary_id":ObjectId(summary_id),"specimen_id": specimen_id},{"$set":{ref:hit_ref_nums}})
        except Exception, e:
            self.bind_object.logger.info("导入%s结果表出错:%s" % (anno, e))
        else:
            self.bind_object.logger.info("导入%s结果表成功!" % anno)