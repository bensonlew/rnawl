# -*- coding: utf-8 -*-
# __author__ = 'zouxuan'
# from mainapp.config.db import get_mongo_client
from bson.objectid import ObjectId
import types
from bson import SON
from mainapp.models.workflow import Workflow
from mainapp.models.mongo.meta import Meta
import random
import json


class FungiGenome(Meta):
    def __init__(self, bind_object=None):
        self._bind_object = bind_object
        super(FungiGenome, self).__init__(self._bind_object)
        self._project_type = "fungigenome"

    def get_projectsn(self, task_id):
        collection = self.db['sg_task']
        result = collection.find_one({"task_id": task_id})
        project_sn = result["project_sn"]
        return project_sn
    def get_genefile_bysample(self, task_id, sample, type="faa", predict = "gene"):  # add by zhujuan 20180409
        if type not in ["fnn", "gff", "faa"]:
            raise Exception("输入type参数必须为fnn, gff, faa!")
        if predict not in ["gene", "trna", "rrna"]:
            raise Exception("输入predict必须为gene, trna, rrna!")
        collection = self.db[predict + '_predict']
        result = collection.find_one({"task_id": task_id})
        result_path = result['file_path']
        sample_gene_path = {}
        for i in sample.split(","):
            file_path = ''.join([result_path[0], i, result_path[1], i, result_path[2], type])
            sample_gene_path[i] = file_path
        return sample_gene_path

    def get_middle_path(self, task_id, predict = 'gene'):
        collection = self.db[predict + '_predict']
        result = collection.find_one({"task_id": task_id})
        result_path = result['file_path']
        #middle_path = 'rerewrweset' + result_path[0].split('rerewrweset')[-1]
        middle_path = '/'.join(result_path[0].rstrip('/').split('/')[-5:])
        #example  rerewrweset/files/m_188/188_5b459756d0ec3/tsg_32207/workflow_results
        return middle_path

    def update_mongo(self,db_name,search, change):
        db = self.db[db_name]
        ret = db.find_one(search)
        if ret:
            db.update({"_id":ret["_id"]},{"$set":change})
            return str(ret["_id"])
        else:
            return 'not find'


