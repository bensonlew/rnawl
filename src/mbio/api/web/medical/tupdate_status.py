# -*- coding: utf-8 -*-
# __author__ = 'moli.zhou'
import json
import datetime
from bson.objectid import ObjectId
from biocluster.core.function import filter_error_info
from ..meta.update_status import UpdateStatus


class TupdateStatus(UpdateStatus):

    def __init__(self, data):
        super(TupdateStatus, self).__init__(data)
        self._client = "client03"
        # self._key = "hM4uZcGs9d"
        # self._url = "http://api.tsanger.com/task/add_file"
        # self.database = self._mongo_client[self.config.MONGODB + '_paternity_test']
        self._binds_id = "5f4324a19b79000093008059"
        self._interface_id = 66
        self._env_name = "offline"
        self._key = "458b97de4c0bb5bf416c8cea208309ed"
        self._url = "http://apicenter.nsg.com/index/in"
        
    def update(self):
        pass

    def update_status(self):
        if "status" in self.data["sync_task_log"]["task"].keys():
            status = self.data["sync_task_log"]["task"]["status"]
        else:
            return
        desc = ''
        for i in self.data['sync_task_log']['log']:
            if 'step_name' not in i:
                desc = i['desc']
        desc = filter_error_info(desc)
        create_time = str(self.data["sync_task_log"]["task"]["created_ts"])

        if not self.update_info:
            return
        for obj_id, collection_name in self.update_info.items():
            obj_id = ObjectId(obj_id)
            collection = self.db[collection_name]
            if status != "start":
                data = {
                    "status": "end" if status == 'finish' else status,
                    "desc": desc,
                    "time": create_time
                }
                collection.find_one_and_update({"_id": obj_id}, {'$set': data}, upsert=True)
 