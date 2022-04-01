# -*- coding: utf-8 -*-
# __author__ = 'HD'
import urllib
import json
import datetime
import re
import sys
from bson.objectid import ObjectId
from biocluster.config import Config as Conf
from biocluster.wpm.log import Log
from biocluster.core.function import CJsonEncoder, filter_error_info
import traceback
import random
import hashlib
import time


class UpdateStatus(Log):
    """
    meta的web api，用于更新sg_status表并向前端发送状态信息和文件上传信息
    一般可web api功能可从此处继承使用，需要重写__init__方法
    """

    def __init__(self, data):
        super(UpdateStatus, self).__init__(data)
        self.config = Conf()
        # self._client = Config().rcf.get("toollab", "client")
        # self._key = Config().rcf.get("toollab", "authkey")
        # self._binds_id = Config().rcf.get("toollab", "binds_id")
        # self._interface_id = Config().rcf.get("toollab", "interface_id")
        # self._env_name = Config().rcf.get("toollab", "env_name")
        webcfg = self.config.get_webauth("tool_lab")
        self._key = webcfg["authkey"]
        self._binds_id = webcfg["binds_id"]
        self._interface_id = webcfg["interface_id"]
        self._env_name = webcfg["env_name"]
        self._url = "http://apicenter.lab.majorbio.com/index/in"
        self._project_db = self.config.get_mongo_client(mtype="project", db_version=self._dbversion)[
            self.config.get_mongo_dbname(mtype="project", db_version=self._dbversion)]

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
            format_task = self.data["sync_task_log"]["task"]
            format_task['task_id'] = self._get_main_task_id()
            content["basis"] = format_task
        if 'files' in self.data['sync_task_log'].keys():
            content['files'] = self.data["sync_task_log"]["files"]
        if 'dirs' in self.data['sync_task_log'].keys():
            content['dirs'] = self.data['sync_task_log']['dirs']
        if 'base_path' in self.data['sync_task_log'].keys():
            content['basis']['base_path'] = self.data['sync_task_log']['base_path']
        if "region" in self.data['sync_task_log'].keys():
            content['basis']['region'] = self.data['sync_task_log']['region']
        if 'log' in self.data['sync_task_log'].keys():
            content['log'] = self.data['sync_task_log']['log']
        content['basis']['project_sn'] = self.get_info_from_maintable('project_sn')
        data['sync_task_log'] = json.dumps(content, cls=CJsonEncoder)
        yield urllib.urlencode(data)

    def get_info_from_maintable(self, key):
        """
        通过查找主表中key对应的value值
        :param key:
        :return:
        """
        resule_data = "None"
        if not self.update_info:
            return
        # noinspection PyBroadException
        try:
            for obj_id, collection_name in self.update_info.items():
                obj_id = ObjectId(obj_id)
                temp_find = self._get_params(collection_name, obj_id)
                if temp_find and key in temp_find.keys():
                    resule_data = temp_find[key]
        except Exception, e:
            exstr = traceback.format_exc()
            print exstr
            print e
            sys.stdout.flush()
        return resule_data

    def update(self):
        # 即时的不发送数据给前端
        if self.get_info_from_maintable('instant'):
            pass
        else:
            self.send()

    def get_sig(self):
        nonce = str(random.randint(1000, 10000))
        timestamp = str(int(time.time()))
        x_list = [self._key, timestamp, nonce]
        x_list.sort()
        sha1 = hashlib.sha1()
        map(sha1.update, x_list)
        sig = sha1.hexdigest()
        signature = {
            "nonce": nonce,
            "timestamp": timestamp,
            "signature": sig,
            "binds_id": self._binds_id,
            "interface_id": self._interface_id,
            "env_name": self._env_name
        }
        return urllib.urlencode(signature)

    def _get_main_task_id(self):
        my_id = re.split('_', self.task_id)
        my_id.pop(-1)
        my_id.pop(-1)
        return "_".join(my_id)

    def _get_params(self, collection_name, main_id):
        try:
            return self.db[collection_name].find_one({"_id": main_id})
        except Exception, e:
            exstr = traceback.format_exc()
            print exstr
            print e
            sys.stdout.flush()

    def update_status(self):
        if "status" in self.data["sync_task_log"]["task"].keys():
            status = self.data["sync_task_log"]["task"]["status"]
        else:
            return
        desc = ''
        error = ""
        for i in self.data['sync_task_log']['log']:
            if 'name' not in i:
                desc = i['desc']
                if "error" in i.keys() and isinstance(i['error'], dict):
                    error = i['error']

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
                    "created_ts": create_time
                }
                collection.update_one({"_id": obj_id}, {'$set': data}, upsert=True)
            sg_status_col = self.db['sg_status']
            if status == "start":
                temp_find = self._get_params(collection_name, obj_id)
                params = ""
                submit_location = ""
                if temp_find and "params" in temp_find.keys():
                    params = temp_find['params']
                    # noinspection PyBroadException
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
                sg_status_col.insert_one(insert_data)
            elif status == "finish":  # 只能有一次finish状态
                insert_data = {
                    "status": 'end',
                    "desc": desc,
                    "time": datetime.datetime.now().strftime("%Y-%m-%d %H:%M:%S")
                }
                sg_status_col.update_one({"table_id": obj_id, "type_name": collection_name},
                                         {'$set': insert_data}, upsert=True)
            else:
                insert_data = {
                    "status": status,
                    "desc": desc,
                    "time": datetime.datetime.now().strftime("%Y-%m-%d %H:%M:%S")
                }
                if error:
                    insert_data["error"] = error
                sg_status_col.update_one({"table_id": obj_id, "type_name": collection_name},
                                         {'$set': insert_data}, upsert=True)
