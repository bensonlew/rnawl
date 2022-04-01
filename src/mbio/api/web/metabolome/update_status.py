# -*- coding: utf-8 -*-
# __author__ = 'guoquan'
# last_modified = shaohua.yuan
# -*- coding: utf-8 -*-
# __author__ = 'qiuping'
from ..meta.update_status import UpdateStatus as Us

import json
from biocluster.core.function import CJsonEncoder, filter_error_info
from bson.objectid import ObjectId
import datetime
import sys
import traceback


class UpdateStatus(Us):

    def __init__(self, data):
        super(UpdateStatus, self).__init__(data)
        self._project_type = 'metabolome'

    # edit zouguanqing 20200304
    def update_status(self):
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
                    if temp_find and "params" in temp_find.keys():
                        params = temp_find['params']
                        try:
                            temp_json = json.loads(temp_find['params'])
                        except:
                            pass
                        else:
                            if isinstance(temp_json, dict) and "submit_location" in temp_json.keys():
                                submit_location = temp_json['submit_location']
                    insert_data = {
                        "table_id": obj_id,
                        "table_name": temp_find["name"] if temp_find and "name" in temp_find.keys() else "",
                        "task_id": self._get_main_task_id(),
                        "run_id": self.task_id,
                        "type_name": collection_name,
                        "params": params,
                        "submit_location": submit_location,
                        "status": "start",
                        "is_new": "new",
                        "desc": desc,
                        "time": datetime.datetime.now().strftime("%Y-%m-%d %H:%M:%S")
                    }

                    # zouguanqing 20200304
                    tmp_json = json.loads(params)
                    add_item = ['metab_set','metabset']
                    for each_i in add_item:
                        if each_i in tmp_json.keys():
                            insert_data['metabset_id'] = tmp_json[each_i]
                    if 'metab_table' in  tmp_json.keys():
                        insert_data['metab_table'] = tmp_json['metab_table']

                    if 'assodata_id' in tmp_json.keys():
                        insert_data['assodata_id'] = tmp_json['assodata_id']


                    sg_status_col.insert_one(insert_data)
            elif status == "finish":  # 只能有一次finish状态
                if not batch_id:
                    insert_data = {
                        "status": 'end',
                        "desc": desc,
                        "time": datetime.datetime.now().strftime("%Y-%m-%d %H:%M:%S")
                    }
                    sg_status_col.update_one({"table_id": obj_id, "type_name": collection_name},
                                             {'$set': insert_data}, upsert=True)
                self.pipe_update(batch_id, collection_name, obj_id, "end", desc)
            else:
                if not batch_id:
                    insert_data = {
                        "status": status,
                        "desc": desc,
                        "time": datetime.datetime.now().strftime("%Y-%m-%d %H:%M:%S")
                    }
                    sg_status_col.update_one({"table_id": obj_id, "type_name": collection_name},
                                             {'$set': insert_data}, upsert=True)
                self.pipe_update(batch_id, collection_name, obj_id, status, desc)
        try:
            self._update_project_db(batch_id)
        except Exception, e:
            exstr = traceback.format_exc()
            print exstr
            print e
            sys.stdout.flush()


