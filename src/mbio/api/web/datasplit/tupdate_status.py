# -*- coding: utf-8 -*-
# __editor__ = "wangzhaoyue"

import re
import json
import urllib
from biocluster.wpm.log import Log
from bson.objectid import ObjectId
from biocluster.core.function import CJsonEncoder, filter_error_info
import traceback
import sys
import random
import time
import hashlib
import copy


class TupdateStatus(Log):

    def __init__(self, data):
        super(TupdateStatus, self).__init__(data)
        self._client = "client03"
        # self._key = "hM4uZcGs9d"
        # self._url = "http://api.tsg.com/task/add_file"
        self._project_type = "datasplit"
        self._binds_id = "5f4324a19b79000093008059"
        self._interface_id = 66
        self._env_name = "offline"
        self._key = "458b97de4c0bb5bf416c8cea208309ed"
        self._url = "http://apicenter.nsg.com/index/in"

    def __del__(self):
        self._mongo_client.close()

    @property
    def post_data(self):
        data = dict()
        content = {
            "basis": {
                "task_id": self._get_main_task_id()
            }
        }
        if "task" in self.data['sync_task_log'].keys():
            format_task = copy.deepcopy(self.data["sync_task_log"]["task"])
            format_task['task_id'] = self._get_main_task_id()
            content["basis"] = format_task
        if 'files' in self.data['sync_task_log'].keys():
            content['files'] = self.data["sync_task_log"]["files"]
        if 'dirs' in self.data['sync_task_log'].keys():
            content['dirs'] = self.data['sync_task_log']['dirs']
        if 'base_path' in self.data['sync_task_log'].keys():
            content['basis']['base_path'] = self.data['sync_task_log']['base_path']
            # noinspection PyBroadException
            try:
                project_sn = self.data['sync_task_log']["base_path"].split('/')[3]
            except Exception, e:
                exstr = traceback.format_exc()
                print exstr
                print e
                sys.stdout.flush()
            else:
                content['basis']["project_sn"] = project_sn
        if "region" in self.data['sync_task_log'].keys():
            content['basis']['region'] = self.data['sync_task_log']['region']
        if 'log' in self.data['sync_task_log'].keys():
            content['log'] = self.data['sync_task_log']['log']
        data['sync_task_log'] = json.dumps(content, cls=CJsonEncoder)
        yield urllib.urlencode(data)

    # @property
    # def post_data(self):
    #     data = dict()
    #     content = {
    #         "task": {
    #             "task_id": self.task_id
    #         }
    #     }
    #     if "files" in self.data["sync_task_log"].keys():
    #         content["files"] = self.data["sync_task_log"]["files"]
    #     if "dirs" in self.data["sync_task_log"].keys():
    #         content["dirs"] = self.data["sync_task_log"]["dirs"]
    #     if "base_path" in self.data["sync_task_log"].keys():
    #         content["base_path"] = self.data["sync_task_log"]["base_path"]
    #     data["sync_task_log"] = json.dumps(content, cls=CJsonEncoder)
    #     return urllib.urlencode(data)

    def get_sig(self):
        """
        增加bind_id等前端需要的验证
        :return:
        """
        nonce = str(random.randint(1000, 10000))
        timestamp = str(int(time.time()))
        x_list = [self._key, timestamp, nonce]
        x_list.sort()
        sha1 = hashlib.sha1()
        map(sha1.update, x_list)
        sig = sha1.hexdigest()
        signature = {
            "client": self._client,
            "nonce": nonce,
            "timestamp": timestamp,
            "signature": sig,
            "binds_id": self._binds_id,
            "interface_id": self._interface_id,
            "env_name": self._env_name
        }
        return urllib.urlencode(signature)

    def update(self):
        pass

    def update_status(self):
        if "status" in self.data["sync_task_log"]["task"].keys():
            status = self.data["sync_task_log"]["task"]["status"]
        else:
            return
        desc = ""
        for i in self.data["sync_task_log"]["log"]:
            if "name" not in i:
                desc = i["desc"]
        end_time = str(self.data["sync_task_log"]["task"]["created_ts"])
        if not self.update_info:
            return
        for obj_id, collection_name in self.update_info.items():
            obj_id = ObjectId(obj_id)
            collection = self.db[collection_name]
            self.logger.info("datasplit:{}".format(collection))
            result = collection.find_one({"_id": obj_id})
            if collection_name == "sg_split":
                split_status = result["split_status"]
                data = {
                    "desc": desc
                }
                if split_status["first_split"] == "start":
                    if status != "start":
                        split_status["first_split"] = "end" if status == "finish" else status
                        data = {
                            "split_status": split_status,
                            "desc": "文库拆分运行结束",
                            "end_ts": end_time
                        }
                        if result["split_model"] == "manual":
                            data["status"] = "end" if status == "finish" else status
                            data["desc"] = desc
                        if status != "finish":
                            data["status"] = status
                            data["desc"] = desc
                elif split_status["second_split"] == "start":
                    if status != "start":
                        split_status["second_split"] = "end" if status == "finish" else status
                        data = {
                            "split_status": split_status,
                            "desc": "样本拆分运行结束",
                            "end_ts": end_time
                        }
                        if result["split_model"] == "manual":
                            data["status"] = "end" if status == "finish" else status
                            data["desc"] = desc
                        if status != "finish":
                            data["status"] = status
                            data["desc"] = desc
                elif split_status["qc"] == "start":
                    if status != "start":
                        split_status["qc"] = "end" if status == "finish" else status
                        data = {
                            "split_status": split_status,
                            "desc": desc,
                            "end_ts": end_time,
                            "status": "end" if status == "finish" else status
                        }
                elif split_status["cpc"] == "start":
                    if status != "start":
                        split_status["cpc"] = "end" if status == "finish" else status
                        data = {
                            "split_status": split_status,
                            "desc": desc,
                            "cpc_ts": end_time,
                            "cpc": "end" if status == "finish" else status
                        }
                print "datasplit_status:data:"
                print split_status
                print data
                collection.find_one_and_update({"_id": obj_id}, {"$set": data}, upsert=True)   # 更新状态表
            else:
                if status != "start":
                    data = {
                        "status": "end" if status == 'finish' else status,
                        "desc": desc,
                        "end_time": end_time
                    }
                    collection.update({"_id": obj_id}, {'$set': data}, upsert=True)
