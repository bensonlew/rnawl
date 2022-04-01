# -*- coding: utf-8 -*-
# __author__ = 'zengjing'
# last_modify: 2016.12.06

from bson.objectid import ObjectId
import datetime
import types
from mainapp.models.mongo.core.base import Base
# from mainapp.config.db import get_mongo_client
# from biocluster.config import Config


class FunctionPredict(Base):
    def __init__(self, bind_object=None):
        super(FunctionPredict, self).__init__(bind_object)
        self._project_type = 'meta'
        # self.client = get_mongo_client()
        # self.db = self.client[Config().MONGODB]

    def add_function_predict(self, name=None, params=None, otu_id=None):
        if otu_id != 0 and not isinstance(otu_id, ObjectId):
            if isinstance(otu_id, types.StringTypes):
                otu_id = ObjectId(otu_id)
            else:
                raise Exception("otu_id!")
        collection = self.db["sg_otu"]
        result = collection.find_one({"_id": otu_id})
        project_sn = result['project_sn']
        task_id = result['task_id']
        insert_data = {
            "project_sn": project_sn,
            "task_id": task_id,
            "name": name if name else "16s_function_predict_" + str(datetime.datetime.now().strftime("%Y%m%d_%H%M%S")),
            "params": params,
            "status": "end",
            "desc": "16s功能预测主表",
            "created_ts": datetime.datetime.now().strftime("%Y-%m-%d %H:%M:%S"),
        }
        collection = self.db["sg_16s_function_prediction"]
        predict_id = collection.insert_one(insert_data).inserted_id
        return predict_id
