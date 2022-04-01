# -*- coding: utf-8 -*-
# __author__ = 'shicaiping'
# __date__ = '20200320

from biocluster.api.database.base import Base, report_check
import os
import datetime
import pymongo
from bson.objectid import ObjectId
from types import StringTypes
from bson.son import SON
from mbio.api.database.ref_rna_v2.api_base import ApiBase


class AnnotationInfo(ApiBase):
    def __init__(self, bind_object=None):
        super(AnnotationInfo, self).__init__(bind_object)
        self._project_type = 'ref_genome_db'

    @report_check
    def add_anno(self, client, anno_stat):
        if client == "client03":
            self._db = self.db
        else:
            tsg_client = pymongo.MongoClient(
                "mongodb://rna:y6a3t5n1w0y7@10.8.0.23/sanger_ref_rna_v2?authMechanism=SCRAM-SHA-1")
            self._db = tsg_client["sanger_ref_genome_db"]
        main_insert_data = {
            'project_sn': self.bind_object.sheet.project_sn,
            'task_id': self.bind_object.sheet.id,
            'created_ts': datetime.datetime.now().strftime('%Y-%m-%d %H:%M:%S'),
            'version': 'v2',
            'status': 'end'
        }

        results = self._db['sg_anno'].find({"task_id": self.bind_object.sheet.id})
        for result in results:
            main_id = result["_id"]
            self.db['sg_anno_detail'].delete_many({"stat_id": main_id})
            self.db['sg_anno'].delete_one(result)

        try:
            main_collection = self._db['sg_anno']
            insert_id = main_collection.insert_one(main_insert_data).inserted_id
            self.add_anno_detail(insert_id, anno_stat)
        except Exception, e:
            self.bind_object.set_error("注释统计信息导入失败, {}".format(e))
        else:
            self.bind_object.logger.info('注释统计信息导入成功')

    @report_check
    def add_anno_detail(self, stat_id, anno_stat):
        if stat_id is None:
            self.bind_object.set_error("主表ID不能为空")
        if not isinstance(stat_id, ObjectId):
            if isinstance(stat_id, StringTypes):
                stat_id = ObjectId(stat_id)
            else:
                self.bind_object.set_error("stat_id必须为ObjectId对象或其对应的字符串!")
        insert_data = list()
        with open(anno_stat, "r") as f:
            head = f.readline()
            for line in f:
                items = line.strip().split("\t")
                data = [
                    ('stat_id', stat_id),
                    ('type', items[0]),
                    ('transcript', int(items[1])),
                    ('gene', int(items[2])),
                    ('transcript_percent', round(float(items[3]), 4)),
                    ('gene_percent', round(float(items[4]), 4)),
                ]
                data = SON(data)
                insert_data.append(data)
        try:
            collection = self._db['sg_anno_detail']
            collection.insert_many(insert_data)
        except Exception, e:
            self.bind_object.set_error("导入注释统计信息sg_anno_detail失败, {}".format(e))
        else:
            self.bind_object.logger.info("导入注释统计信息sg_anno_detail成功")