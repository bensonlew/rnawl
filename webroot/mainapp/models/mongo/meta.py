# -*- coding: utf-8 -*-
# __author__ = 'guoquan'
# last_modified = "guhaidong"
from .core.base import Base
# from mainapp.config.db import get_mongo_client
from bson.objectid import ObjectId
import types
# from biocluster.config import Config
from bson import SON
import json
import re
import datetime


class Meta(Base):
    def __init__(self, bind_object=None):
        self._bind_object = bind_object
        super(Meta, self).__init__(self._bind_object)
        self._project_type = "meta"
        # self.client = get_mongo_client()
        '''
        if not db:
            self._db_name = "meta"
            # self.db = self.client[Config().MONGODB]
        else:
            self._db_name = db
            # self.db = self.client[db]
        '''

    '''
    def __del__(self):
        self.client.close()
    '''

    def get_otu_table_info(self, otu_id):

        if isinstance(otu_id, types.StringTypes):
            otu_id = ObjectId(otu_id)
        elif isinstance(otu_id, ObjectId):
            otu_id = otu_id
        else:
            raise Exception("输入otu_id参数必须为字符串或者ObjectId类型!")
        collection = self.db['sg_otu']
        result = collection.find_one({"_id": otu_id})
        return result

    def get_task_info(self, task_id):
        collection = self.db['sg_task']
        result = collection.find_one({"task_id": task_id})
        return result

    def insert_main_table(self, collection, data):
        return self.db[collection].insert_one(SON(data)).inserted_id

    def insert_none_table(self, collection):
        return self.db[collection].insert_one({}).inserted_id

    def insert_main_table_new(self, collection, obj_id, data):
        # return self.db[collection].find_one_and_update({"_id": ObjectId(obj_id)}, {'$set': data}, upsert=True)
        return self.db[collection].update({"_id": ObjectId(obj_id)}, {'$set': data}, upsert=True)

    def update_status_failed(self, collection, doc_id):
        """
        改特定_id主表的status状态从start为failed，主要用于特殊投递任务失败

        params collection: 主表collection名称
        params doc_id: 主表_id
        """
        self.db[collection].update_one({'_id': ObjectId(doc_id), "status": "start"}, {"$set": {'status': 'failed'}})

    def update_workflow_id(self, collection, main_id, workflow_id):
        """
        """
        self.db.workflowid2analysisid.insert_one({'workflow_id': workflow_id,
                                                  "main_id": ObjectId(main_id),
                                                  'collection': collection,
                                                  "created_ts": datetime.datetime.now().strftime("%Y-%m-%d %H:%M:%S")})

    def sampleIdToName(self, sampleIds):
        """
        将一个用逗号隔开的样本ID的集合转换成样本名，返回一个用逗号隔开的样本名的集合
        """
        myIds = re.split("\s*,\s*", sampleIds)
        collection = self.db["sg_specimen"]
        mySampleNames = list()
        for id_ in myIds:
            if id_ == "":
                raise Exception("存在空的sample_id")
            if not isinstance(id_, ObjectId):
                if isinstance(id_, types.StringTypes):
                    id_ = ObjectId(id_)
                else:
                    raise Exception("样本id必须为ObjectId对象或者其对应的字符串！")
            result = collection.find_one({"_id": id_})
            if not result:
                raise Exception("无法根据传入的_id:{}在sg_speciem表里找到相应的记录".format(str(id_)))
            mySampleNames.append(result["specimen_name"])
        mySamples = ",".join(mySampleNames)
        return mySamples

    def get_main_info(self, main_id, collection_name):
        if isinstance(main_id, types.StringTypes):
            main_id = ObjectId(main_id)
        elif isinstance(main_id, ObjectId):
            main_id = main_id
        else:
            raise Exception("输入main_id参数必须为字符串或者ObjectId类型!")
        collection = self.db[collection_name]
        express_info = collection.find_one({'_id': main_id})
        return express_info

    def get_mongo_common(self, db_name, search):
        collection = self.db[db_name]
        result = collection.find_one(search)
        if result:
            return result
        else:
            return False

    def get_new_env_detail(self, env_detail, env_id):
        env_info = self.get_mongo_common('sg_env', ObjectId(env_id))
        env_detail = json.loads(env_detail)
        if 'env_true_name' in env_info:
            true_env = env_info['env_true_name']
        else:
            true_env = None
        if true_env:
            env_names = env_info['env_names'].split(',')
            env_name_dict = dict(zip(true_env, env_names))
            for k in env_detail:
                for i in range(len(env_detail[k])):
                    n = env_detail[k][i]
                    env_detail[k][i] = env_name_dict[n]
        return json.dumps(env_detail)

    def get_new_env_units(self, units, env_id):
        env_info = self.get_mongo_common('sg_env', ObjectId(env_id))
        try:
            units = json.loads(units)
        except ValueError as e:  # 兼容多样性202010阶段性升级前的传参方式
            return units
        if 'env_true_name' in env_info:
            true_env = env_info['env_true_name']
        else:
            true_env = None
        if true_env:
            env_names = env_info['env_names'].split(',')
            env_name_dict = dict(zip(true_env, env_names))
            for i in range(len(units)):
                n = units[i]
                units[i] = env_name_dict[n]
        return ','.join(units)

    def get_bugbase_info(self, bugbase_id):

        if isinstance(bugbase_id, types.StringTypes):
            bugbase_id = ObjectId(bugbase_id)
        elif isinstance(bugbase_id, ObjectId):
            bugbase_id = bugbase_id
        else:
            raise Exception("输入bugbase_id参数必须为字符串或者ObjectId类型!")
        collection = self.db['sg_bugbase']
        result = collection.find_one({"_id": bugbase_id})
        return result

    def get_faprotax_info(self, faprotax_id):

        if isinstance(faprotax_id, types.StringTypes):
            faprotax_id = ObjectId(faprotax_id)
        elif isinstance(faprotax_id, ObjectId):
            faprotax_id = faprotax_id
        else:
            raise Exception("输入faprotax_id参数必须为字符串或者ObjectId类型!")
        collection = self.db['sg_faprotax']
        result = collection.find_one({"_id": faprotax_id})
        return result