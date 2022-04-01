# -*- coding: utf-8 -*-
# __author__ = 'guoquan'

import traceback
import json
from biocluster.wpm.client import worker_client, wait
from biocluster.config import Config
from bson.objectid import ObjectId
import re
import web


class Basic(object):
    def __init__(self, data, instant=False):
        self.dbversion = self._get_dbversion()
        self._instant = instant
        self._json = data
        self._id = data["id"]
        self._return_msg = None
        self._project_type = self.__check_project_type()
        if self._project_type is not None and self._project_type not in ["paternity_test"]:
            self.client = Config().get_mongo_client(mtype=self._project_type, db_version=self.dbversion)
            self.db = self.client[Config().get_mongo_dbname(self._project_type, db_version=self.dbversion)]
        elif "db_type" in data.keys() and data['db_type'] is not None:
            self._project_type = data['db_type']
            self.client = Config().get_mongo_client(mtype=self._project_type, db_version=self.dbversion)
            self.db = self.client[Config().get_mongo_dbname(self._project_type, db_version=self.dbversion)]
        else:
            if "main_table_data" in self._json['options'].keys():
                raise Exception("请输入数据库类型:{},并请保持与main_conf中的名字一致！".format(self._project_type))

    @property
    def id(self):
        """
        获取运行任务的ID

        :return:
        """
        return self._id

    @property
    def instant(self):
        """
        任务是否是即时计算

        :return: bool
        """
        return self._instant

    @property
    def return_msg(self):
        """
        获取运行任务的返回值

        :return:
        """
        return self._return_msg

    def __check_project_type(self):
        module_name = self.__module__
        mlist = module_name.split(".")
        mlist.pop()
        m = re.match(r'mainapp.controllers.(?:instant|submit).([\w_]+)', ".".join(mlist))
        if m:
            return m.group(1)
        else:
            return None

    def run(self):
        data = web.input()
        dbversion = self._get_dbversion()
        dydb = self._get_dydb()
        contract_id = self._get_contract_id()
        dydb_contract_ids = self._get_dydb_contract_ids()
        self._json['creator'] = "webSubmit"
        if self._project_type is not None:
            self._json['creator'] = "webSubmitBy"+self._project_type
        # print "dbversion:{}".format(dbversion)
        if "main_table_data" in self._json['options'].keys():
            main_table_data = self._json['options']['main_table_data']
        else:
            main_table_data = None
        if 'update_info' in self._json['options'].keys():
            update_info = json.loads(self._json["options"]['update_info'])
            if "update_info" not in self._json.keys():  # by zengjing 20211220 解决data.json里的update_info为空问题
                self._json["update_info"] = self._json["options"]['update_info']
        else:
            update_info = None
        if "main_table_data" in self._json['options'].keys():
            del self._json['options']['main_table_data']
        if dbversion != 0:
            self._json['DBVersion'] = dbversion
        has_id, test_id = self._get_task_test_id()
        if has_id:
            self._json['user_tasks_test_id'] = test_id
        if hasattr(data, "wfm_port"):
            self._json["wfm_port"] = data.wfm_port
        if "stage_id" in self._json.keys():
            try:
                self._json["stage_id"] = str(self._json["stage_id"])
            except:
                self._json["stage_id"] = "0"
        if "project_sn" not in self._json.keys():
            self._json["project_sn"]=self._get_projectsn()
        try:
            if isinstance(self._json["project_sn"], list):
                self._json["project_sn"] = str(self._json["project_sn"][0])
        except:
            pass
        if dydb !=0:
            self._json['dydb'] = dydb
            if contract_id == "":
                return {"success": False, "info": "dydb为1时必须设置contract_id"}
            self._json['contract_id'] = contract_id
        if dydb_contract_ids:
            self._json['dydb_contract_ids'] = dydb_contract_ids
        if update_info:
            for i in update_info:
                if i == "batch_id":
                    continue
                # self.insert_main_table_new(update_info[i], i, {"origin_task_id": self._json['id']})
                if main_table_data:
                    self.insert_main_table_new(update_info[i], i, main_table_data)
        try:
            worker = worker_client()
            info = worker.add_task(self._json)
            if "success" in info.keys() and info["success"]:
                # if update_info:
                #     for i in update_info:
                #         if i == "batch_id":
                #             continue
                #         # self.insert_main_table_new(update_info[i], i, {"origin_task_id": self._json['id']})
                #         if main_table_data:
                #             self.insert_main_table_new(update_info[i], i, main_table_data)
                return info
            else:
                return {"success": False, "info": "任务提交失败 %s" % (info["info"])}
        except Exception, e:
            exstr = traceback.format_exc()
            print "ERROR:", exstr
            raise Exception("任务提交失败：%s, %s" % (str(e), str(exstr)))

    def insert_main_table_new(self, collection, obj_id, data):
        return self.db[collection].update({"_id": ObjectId(obj_id)}, {'$set': data}, upsert=True)

    @staticmethod
    def _get_dbversion():
        data = web.input()
        dbversion = 0
        if hasattr(data, "db_version"):
            # print "data has db_version"
            dbversion = int(data.db_version)
        return dbversion

    @staticmethod
    def _get_dydb():
        data = web.input()
        dydb = 0
        if hasattr(data, "dydb"):
            # print "data has db_version"
            dydb = int(data.dydb)
        return dydb

    @staticmethod
    def _get_projectsn():
        data = web.input()
        project_sn = ""
        if hasattr(data, "project_sn"):
            # print "data has db_version"
            project_sn = data.project_sn
        return project_sn

    @staticmethod
    def _get_contract_id():
        data = web.input()
        contract_id = ""
        if hasattr(data, "contract_id"):
            contract_id = data.contract_id
        return contract_id

    @staticmethod
    def _get_dydb_contract_ids():
        data = web.input()
        dydb_contract_ids = {}
        if hasattr(data, "dydb_contract_ids"):
            # dydb_contract_ids = data.dydb_contract_ids
            try:
                dydb_contract_ids = json.loads(data.dydb_contract_ids)
            except Exception as e:
                print(e)
                dydb_contract_ids = data.dydb_contract_ids
        return dydb_contract_ids

    @staticmethod
    def _get_task_test_id():
        data = web.input()
        if hasattr(data, "user_tasks_test_id"):
            try:
                uttid = int(data.user_tasks_test_id)
            except:
                uttid = 0
            return True, uttid
        else:
            return False, 0
