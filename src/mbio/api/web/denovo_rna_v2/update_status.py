# -*- coding: utf-8 -*-
import json
import datetime
import sys
from bson.objectid import ObjectId
from biocluster.core.function import CJsonEncoder, filter_error_info
import traceback
from ..meta.update_status import UpdateStatus as Us


class UpdateStatus(Us):

    def __init__(self, data):
        super(UpdateStatus, self).__init__(data)
        self._project_type = "denovo_rna_v2"

    def update_status(self):
        super(UpdateStatus, self).update_status()
        if "status" in self.data["sync_task_log"]["task"].keys():
            status = self.data["sync_task_log"]["task"]["status"]
        else:
            return
        desc = ''
        for i in self.data['sync_task_log']['log']:
            if 'name' not in i:
                desc = i['desc']
        desc = filter_error_info(desc)
        create_time = str(self.data["sync_task_log"]["task"]["created_ts"])
        if not self.update_info:
            return
        batch_id = self.update_info.pop("batch_id") if 'batch_id' in self.update_info else None
        print self.update_info
        for obj_id, collection_name in self.update_info.items():
            obj_id = ObjectId(obj_id)
            collection = self.db[collection_name]
            if status != "start":
                data = {
                    "status": "end" if status == 'finish' else status,
                    "desc": desc,
                    "created_ts": create_time
                }
                collection.update_one({"_id": obj_id}, {'$set': data}, upsert=True)
            sg_status_col = self.db['sg_status']
            if status == "start":
                if not batch_id:
                    temp_find = self._get_params(collection_name, obj_id)
                    params = ""
                    submit_location = ""
                    group_id = ""
                    control_id = ""
                    geneset_id = ""
                    if temp_find and "params" in temp_find.keys():
                        params = temp_find['params']
                        try:
                            temp_json = json.loads(temp_find['params'])
                        except:
                            pass
                        else:
                            if isinstance(temp_json, dict) and "submit_location" in temp_json.keys():
                                submit_location = temp_json['submit_location']
                            if isinstance(temp_json, dict) and "group_id" in temp_json.keys():
                                group_id = temp_json['group_id']
                            if isinstance(temp_json, dict) and "control_id" in temp_json.keys():
                                control_id = temp_json['control_id']
                            if isinstance(temp_json, dict) and "geneset_id" in temp_json.keys():
                                geneset_id = temp_json['geneset_id']
                    insert_data = {
                        # "table_id": obj_id,
                        # "table_name": temp_find["name"] if temp_find and "name" in temp_find.keys() else "",
                        # "task_id": self._get_main_task_id(),
                        # "type_name": collection_name,
                        # "params": params,
                        # "submit_location": submit_location,
                        # "status": "start",
                        # "is_new": "new",
                        "group_id": group_id,
                        "control_id": control_id,
                        "geneset_id": geneset_id,
                        # "desc": desc,
                        # "time": datetime.datetime.now().strftime("%Y-%m-%d %H:%M:%S")
                    }
                    sg_status_col.update_one({"table_id": obj_id, "type_name": collection_name},
                                             {'$set': insert_data}, upsert=True)