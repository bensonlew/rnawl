# -*- coding: utf-8 -*-
# __author__ = 'xuting'

from mainapp.models.mongo.core.base import Base
import datetime

# from mainapp.config.db import get_mongo_client
# from biocluster.config import Config


class SampleExtract(Base):
    def __init__(self, bind_object=None):
        super(SampleExtract, self).__init__(bind_object)
        self._project_type = 'meta'
        # self.client = get_mongo_client()
        # self.db = self.client[Config().MONGODB]

    def add_sg_seq_sample(self, task_id, file_path, params, query_id):
        insert_data = {
            "task_id": task_id,
            "file_path": file_path,
            "params": params,
            "query_id": query_id
        }
        collection = self.db["sg_seq_sample"]
        try:
            results = collection.find({"task_id": task_id, "query_id": query_id}).sort("_id",-1)[0] ##add by qingchen.zhang @20191219 目的是为了防止前面的判断出错，导致sg_seq_sample多插入一条记录
            if not results:
                inserted_id = collection.insert_one(insert_data).inserted_id
            else:
                inserted_id = results["_id"]
        except:
            inserted_id = collection.insert_one(insert_data).inserted_id
        return inserted_id

    def add_sample_check(self, task_id, file_path, params, query_id):
        """
        多样性导入样本检测主表
        :param task_id: 任务id
        :param file_path: info_path
        :param params: params
        :param query_id: 文件或者文件夹的id
        :return:
        """
        insert_data = {
            "task_id": task_id,
            "file_path": file_path,
            "params": params,
            "query_id": query_id,
            'status': 'start',
            'desc': 'sample_check主表',
            'created_ts': datetime.datetime.now().strftime("%Y-%m-%d %H:%M:%S"),
        }
        collection = self.db["sg_sample_check"]
        inserted_id = collection.insert_one(insert_data).inserted_id
        return inserted_id

    def add_main_table(self, task_id, from_task, params, query_id):
        """
        样本检测，copy情况插入主表
        :param task_id:
        :param from_task:
        :param params:
        :param query_id:
        :return:
        """
        insert_data = {
            "task_id": task_id,
            "from_task": from_task,
            "params": params,
            "query_id": query_id,
            "status": "start",
            'created_ts': datetime.datetime.now().strftime("%Y-%m-%d %H:%M:%S")
        }
        collection = self.db["sg_sample_check"]
        inserted_id = collection.insert_one(insert_data).inserted_id
        return inserted_id
