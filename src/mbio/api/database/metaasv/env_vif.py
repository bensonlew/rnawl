# -*- coding: utf-8 -*-
# __author__ = 'qingchen.zhang'


from biocluster.api.database.base import Base, report_check
import json
import datetime
from bson.objectid import ObjectId


class EnvVif(Base):
    def __init__(self, bind_object):
        super(EnvVif, self).__init__(bind_object)
        self._project_type = "metaasv"

    @report_check
    def add_env_vif(self, anno_type, name=None, params=None):
        task_id = self.bind_object.sheet.id
        insert_data = {
            "project_sn": self.bind_object.sheet.project_sn,
            "task_id": task_id,
            "status": "end",
            "desc": "vif_Origin",
            "params": json.dumps(params, sort_keys=True, separators=(',', ':')),
            "created_ts": datetime.datetime.now().strftime("%Y-%m-%d %H:%M:%S")
        }
        collection = self.db["vif"]
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
            if table_type == "before":
                before_vif_list = envs_list
            else:
                after_vif_list = envs_list
            vif_value = f.readline().strip().split("\t")
            data = {
                    "name": "VIF",
                    "type": table_type,
                    "vif_id": ObjectId(vif_id),
                   }
            for n, e in enumerate(envs):
                data[e] = vif_value[n]
            data_list.append(data)
        try:
            collection = self.db["vif_detail"]
            collection.insert_many(data_list)
            self.bind_object.logger.info("导入detail表成功")
            main_collection = self.db["vif"]
            settled_params = {"software" : "R-3.3.1 (vegan)"}
            settled_params_json = json.dumps(settled_params, sort_keys=True, separators=(',', ':'))
            if table_type == "before":
                before_table_data = {
                "table_data": envs,
                            "condition": {"type": "before"}}
                table_data_json = json.dumps(before_table_data, sort_keys=True, separators=(',', ':'))
                main_collection.update({"_id": ObjectId(vif_id)}, {"$set": {"before_env_list": before_vif_list,
                                                                            "main_id": ObjectId(vif_id),
                                                                            "settled_params": settled_params_json,
                                                                            "before_table_data": table_data_json}})
            else:
                after_table_data = {
                "table_data": envs,
                            "condition": {"type": "after"}}
                table_data_json = json.dumps(after_table_data, sort_keys=True, separators=(',', ':'))
                main_collection.update({"_id": ObjectId(vif_id)}, {"$set": {"after_env_list": after_vif_list,
                                                                            "main_id": ObjectId(vif_id),
                                                                            "settled_params": settled_params_json,
                                                                            "after_table_data": table_data_json}})
            self.bind_object.logger.info("更新主表成功")
        except Exception, e:
            self.bind_object.logger.error("导入env_vif结果数据出错:%s" % e)
            self.bind_object.set_error("导入env_vif数据出错")
        else:
            self.bind_object.logger.info("导入env_vif结果数据成功")
