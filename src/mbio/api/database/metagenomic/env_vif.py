# -*- coding: utf-8 -*-
# __author__ = 'zhujuan'


from biocluster.api.database.base import Base, report_check
import json
import datetime
from bson.objectid import ObjectId
#from biocluster.config import Config


class EnvVif(Base):
    def __init__(self, bind_object):
        super(EnvVif, self).__init__(bind_object)
        self._project_type = "metagenomic"
        #self._db_name = Config().MONGODB + '_metagenomic'

    @report_check
    def add_env_vif(self, anno_type, name=None, params=None):
        task_id = self.bind_object.sheet.id
        insert_data = {
            "project_sn": self.bind_object.sheet.project_sn,
            "task_id": task_id,
            "anno_type": anno_type,
            #"env_id": env_id,
            #"geneset_id": geneset_id,
            #"anno_id": anno_id,
            "name": name,
            #"level_id": level_id,
            "status": "end",
            #'group_id': group_id,
            "desc": "vif方差膨胀因子分析",
            "params": json.dumps(params, sort_keys=True, separators=(',', ':')),
            "created_ts": datetime.datetime.now().strftime("%Y-%m-%d %H:%M:%S")
        }
        collection = self.db["env_vif"]
        inserted_id = collection.insert_one(insert_data).inserted_id
        return inserted_id

    @report_check
    def add_env_vif_detail(self, file_path, table_type, vif_id):
        data_list = []
        before_vif_list = ""
        after_vif_list = ""
        with open(file_path, "r") as f:
            envs = f.readline().strip().split("\t")
            envs_list = ','.join(envs)
            if table_type == "before_vif":
                before_vif_list = envs_list
            else:
                after_vif_list = envs_list
            vif_value = f.readline().strip().split("\t")
            data = {
                    "table_type": table_type,
                    "vif_id": ObjectId(vif_id),
                   }
            for n, e in enumerate(envs):
                data[e] = vif_value[n]
            data_list.append(data)
        try:
            collection = self.db["env_vif_detail"]
            collection.insert_many(data_list)
            self.bind_object.logger.info("开始刷新主表写入")
            main_collection = self.db["env_vif"]
            if table_type == "before_vif":
                main_collection.update({"_id": ObjectId(vif_id)}, {"$set": {"before_env_list": before_vif_list}})
            else:
                main_collection.update({"_id": ObjectId(vif_id)}, {"$set": {"after_env_list": after_vif_list}})
        except Exception, e:
            self.bind_object.logger.info("导入env_vif结果数据出错:%s" % e)
            self.bind_object.set_error("导入env_vif结果数据出错", code="52801001")
        else:
            self.bind_object.logger.info("导入env_vif结果数据成功")
