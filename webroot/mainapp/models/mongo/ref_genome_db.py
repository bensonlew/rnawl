# -*- coding: utf-8 -*-
from __future__ import print_function
from bson.objectid import ObjectId
import pymongo
from bson import SON
from mainapp.models.workflow import Workflow
from .core.base import Base
import random
import datetime
__author__ = 'scp'


class RefGenomeDb(Base):
    def __init__(self, bind_object=None):
        self._bind_object = bind_object
        super(RefGenomeDb, self).__init__(bind_object=self._bind_object)
        self._project_type = 'ref_genome_db'

    def delete_genome_db_table(self, client, genome_id):
        if client == "client03":
            self._db = self.db
        else:
            tsg_client = pymongo.MongoClient(
                "mongodb://rna:y6a3t5n1w0y7@10.8.0.23/sanger_ref_rna_v2?authMechanism=SCRAM-SHA-1")
            self._db = tsg_client["sanger_ref_rna_v2"]
        table = self._db["sg_genome_db"]
        result = table.find_one({"genome_id": genome_id})
        if result:
            table.remove({"genome_id": genome_id})

    def get_task_info(self, client, task_id):
        """
        根据task_id到sg_task获得相关记录信息
        :param task_id:
        :return:
        """
        if client == "client03":
            self._db = self.db
        else:
            tsg_client = pymongo.MongoClient(
                "mongodb://rna:y6a3t5n1w0y7@10.8.0.23/sanger_ref_rna_v2?authMechanism=SCRAM-SHA-1")
            self._db = tsg_client["sanger_ref_genome_db"]
        collection = self._db['sg_task']
        result = collection.find({"task_id": task_id}).sort([("_id", -1)])[0]
        return result

    def change_task_status(self, client, task_id):
        """
        根据task_id到sg_task获得相关记录信息
        更改task_status
        """
        if client == "client03":
            self._db = self.db
        else:
            tsg_client = pymongo.MongoClient(
                "mongodb://rna:y6a3t5n1w0y7@10.8.0.23/sanger_ref_rna_v2?authMechanism=SCRAM-SHA-1")
            self._db = tsg_client["sanger_ref_genome_db"]
        collection = self._db['sg_task']
        genome_id = collection.find({"task_id": task_id}).sort([("_id", -1)])[0]["genome_id"]
        self._db['sg_task'].update({'task_id': task_id, 'genome_id': genome_id}, {'$set': {'task_status': "updating..."}}, upsert=True)
        return True

    def rollback_status(self, task_id):
        """
        根据task_id到sg_task获得相关记录信息
        更改task_status
        """
        tsg_client = pymongo.MongoClient(
            "mongodb://rna:y6a3t5n1w0y7@10.8.0.23/sanger_ref_rna_v2?authMechanism=SCRAM-SHA-1")
        self._db = tsg_client["sanger_ref_genome_db"]
        collection = self._db['sg_task']
        genome_id = collection.find({"task_id": task_id}).sort([("_id", -1)])[0]["genome_id"]
        result_dir = collection.find({"task_id": task_id}).sort([("_id", -1)])[0]["result_dir"]
        if 'isanger' in result_dir:
            status = 'isanger_finish'
        elif 'sanger' in result_dir:
            status = 'sanger_finish'
        else:
            status = 'tsg_finish'
        self._db['sg_task'].update({'task_id': task_id, 'genome_id': genome_id}, {'$set': {'task_status': status}}, upsert=True)
        return True

    def insert_main_table(self, client, collection_name, data):
        """
        新建主表或往主表插入一条记录，并返回该条记录的_id
        :param collection_name:
        :param data:
        :return: 插入记录的id
        """
        if client == "client03":
            self._db = self.db
        else:
            tsg_client = pymongo.MongoClient(
                "mongodb://rna:y6a3t5n1w0y7@10.8.0.23/sanger_ref_rna_v2?authMechanism=SCRAM-SHA-1")
            self._db = tsg_client["sanger_ref_genome_db"]
        main_id = self._db[collection_name].insert_one(SON(data)).inserted_id
        conn = self._db[collection_name]
        task_id = conn.find_one({'_id': main_id})['task_id']
        conn.update({'_id': main_id, "task_id": task_id}, {"$set": {'main_id': main_id}}, upsert=True)
        return main_id

    def update_task_info(self, client, task_id):
        if client == "client03":
            self._db = self.db
        else:
            tsg_client = pymongo.MongoClient(
                "mongodb://rna:y6a3t5n1w0y7@10.8.0.23/sanger_ref_rna_v2?authMechanism=SCRAM-SHA-1")
            self._db = tsg_client["sanger_ref_genome_db"]
        collection = self._db['sg_task']
        result = collection.find_one({"task_id": task_id})
        result["task_status"] = "tsg_delete"
        collection.update({"task_id": task_id}, {"$set": result})
        return True

    def get_new_id(self, task_id, otu_id=None):
        """
        根据旧的ID生成新的workflowID，固定为旧的后面用“_”，添加两次随机数或者一次otu_id一次随机数
        """
        if otu_id:
            new_id = "{}_{}_{}".format(task_id, otu_id[-4:], random.randint(1, 10000))
        else:
            # id_ = '%f' % time.time()
            # ids = str(id_).strip().split(".")
            # new_id = "{}_{}_{}".format(task_id, ids[0][5:], ids[1])  #改成时间来命名workflow id
            new_id = "{}_{}_{}".format(task_id, datetime.datetime.now().strftime("%m%d%H%M%S%f"),
                                       random.randint(1, 10000))
        workflow_module = Workflow()
        workflow_data = workflow_module.get_by_workflow_id(new_id)
        if len(workflow_data) > 0:
            return self.get_new_id(task_id, otu_id)
        return new_id


if __name__ == "__main__":
    print('no test Now')
