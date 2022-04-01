# -*- coding: utf-8 -*-
# __author__ = 'gaohao'

from biocluster.api.database.base import Base, report_check
import datetime
import json
from bson.son import SON
import pandas as pd
from bson.objectid import ObjectId


class AnnoSummary(Base):
    def __init__(self, bind_object):
        super(AnnoSummary, self).__init__(bind_object)
        self._project_type = "metagbin"

    @report_check
    def add_anno_summary(self, params=None, name=None, task_id=None):
        if task_id is None:
            task_id = self.bind_object.sheet.id
        project_sn = self.bind_object.sheet.project_sn
        created_ts = datetime.datetime.now().strftime('%Y-%m-%d %H:%M:%S')
        insert_data = {
            "project_sn": project_sn,
            "task_id": task_id,
            "status": "end",
            "desc": "注释汇总",
            "params": json.dumps(params, sort_keys=True, separators=(',', ':')),
            "created_ts": created_ts,
            "name": name if name else 'AnnoSummary_Origin',
        }
        collection = self.db["anno_summary"]
        inserted_id = collection.insert_one(insert_data).inserted_id
        collection.update({'_id': inserted_id},{'$set':{'main_id': inserted_id}})
        return inserted_id

    @report_check
    def add_anno_summary_detail(self, inserted_id, genome_id, anno, database_version='new'):
        data_list = []
        ann = pd.read_table(anno, sep='\t', header=0)
        if len(ann) < 1:
            return
        for i in range(len(ann)):
            data = {
                "summary_id": ObjectId(inserted_id),
                "genome_id": genome_id,
                "gene_id":  ann["Gene ID"][i],
                "location": ann["Location"][i],
                "strand": ann["Strand"][i],
                "start": ann["Start"][i],
                "end": ann["End"][i],
                "len": ann["Gene Length(bp)"][i],
                "gene_des": ann["NR Description"][i],
                "gene_name": ann["Gene Name"][i],
                "cog": ann["COG ID"][i],
                "cog_type": ann["COG Type"][i],
                "ko": ann["KO ID"][i],
                "ko_des": ann["KO Description"][i],
                "cazy_class": ann["Family"][i],
                "cazy_family": ann["Class"][i],
                "cazy_cl_des": ann["Class description"][i],
            }
            if "ARO_name" in list(ann.columns):
                data["aro_name"] = ann["ARO_name"][i]
                data["aro_des"] = ann["ARO_description"][i]
                if database_version in ['new']:
                    data["aro_drug"] = ann["Drug_class"][i]
                    data["aro_resis"] = ann["Resistance_mechanism"][i]
                else:
                    data["aro_category"] = ann["ARO_category"][i]
            data_son = SON(data)
            data_list.append(data_son)
        try:
            collection = self.db["anno_summary_detail"]
            collection.insert_many(data_list)
        except Exception, e:
            self.bind_object.logger.info("导入%s结果表出错:%s" % (anno, e))
        else:
            self.bind_object.logger.info("导入%s结果表成功!" % anno)