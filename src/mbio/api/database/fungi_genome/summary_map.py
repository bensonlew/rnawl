# -*- coding: utf-8 -*-
# __author__ = 'zouxuan'
# last_modify:20180409
from biocluster.api.database.base import Base, report_check
import os
import datetime
import types
import json
from bson.son import SON
from bson.objectid import ObjectId
from biocluster.config import Config


class SummaryMap(Base):
    def __init__(self, bind_object):
        super(SummaryMap, self).__init__(bind_object)
        self._project_type = "fungigenome"
        self.client = Config().get_mongo_client(mtype="metagenomic", ref=True)
        self.mongodb = self.client[Config().get_mongo_dbname("metagenomic", ref=True)]

    @report_check
    def add_map(self, task_id=None):
        kegg_main = self.db['anno_kegg']
        kegg_id = kegg_main.find_one({"task_id":task_id})['_id']
        cazy_main = self.db['anno_cazy']
        cazy_id = cazy_main.find_one({"task_id":task_id})['_id']
        summary_main = self.db['anno_summary']
        sumamary_id = summary_main.find_one({"task_id":task_id})['_id']
        summary_collection = self.db['anno_summary_detail']
        kegg_collection = self.db['anno_kegg_detail']
        cazy_collection = self.db['anno_cazy_detail']
        for i in kegg_collection.find({"kegg_id": kegg_id}):
            gene_id = i['gene_id']
            specimen_id = i['specimen_id']
            map = i['pathway']
            summary_collection.update_one({'summary_id': sumamary_id, 'gene_id': gene_id, 'specimen_id': specimen_id},
                                          {'$set': {'map_id': map}},
                                          upsert=False)
        for i in cazy_collection.find({"cazy_id": cazy_id}):
            gene_id = i['gene_id']
            specimen_id = i['specimen_id']
            class_name = i['class']
            summary_collection.update_one({'summary_id': sumamary_id, 'gene_id': gene_id, 'specimen_id': specimen_id},
                                          {'$set': {'cazy_class': class_name}},
                                          upsert=False)