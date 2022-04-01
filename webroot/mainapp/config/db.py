# -*- coding: utf-8 -*-
# __author__ = 'guoquan'
import web
import os
from biocluster.config import Config as Conf
from biocluster.core.singleton import singleton
import re


@singleton
class Config(object):
    def __init__(self):
        self.config = Conf()
        self.CON_FILE = "./clientkey.txt"
        self.CON_FILE_DICT = {
            "client01": {
                "key": "1ZYw71APsQ",
                "ipmask": "172.16.3.0/24;192.168.10.0/24;127.0.0.1",
                "timelimit": "600",
                "max_workflow": "500"
            },
            "client02": {
                "key": "8A2q9C35Ts",
                "ipmask": "172.16.3.0/24;192.168.10.0/24",
                "timelimit": "",
                "max_workflow": "5"
            },
            "client03": {
                "key": "hM4uZcGs9d",
                "ipmask": "172.16.3.0/24;192.168.10.0/24;127.0.0.1",
                "timelimit": "600",
                "max_workflow": "500"
            },
            "batch": {
                "key": "UD20aMAdga",
                "ipmask": "",
                "timelimit": "",
                "max_workflow": ""
            },
            "rgw": {
                "key": "d2da2f2ca35Gea",
                "ipmask": "",
                "timelimit": "",
                "max_workflow": ""
            },
            "n_client03": {
                "key": "458b97de4c0bb5bf416c8cea208309ed",
                "ipmask": "",
                "timelimit": "",
                "max_workflow": ""
            },
            "n_client01": {
                "key": "e0bb04dd111155b2e6bc6db26d0e1fef",
                "ipmask": "",
                "timelimit": "",
                "max_workflow": ""
            }
        }        
        self.init_config()
        self.RGW_ENABLE = self.config.RGW_ENABLE   # 这里默认都启动对象存储

    def init_config(self):
        """
        检测配置文件clientkey.txt是否存在，如果存在的话就进行初始化
        :return:
        """
        initcfg = {}
        if os.path.exists(self.CON_FILE):
            with open(self.CON_FILE, 'r') as r:
                data = r.readlines()[1:]
                for line in data:
                    tep = line.strip().split("\t")
                    initcfg[tep[1]] = {
                        "key": tep[2],
                        "ipmask": tep[3],
                        "timelimit": tep[4],
                        "max_workflow": tep[5]
                    }
                self.CON_FILE_DICT.update(initcfg)

    def get_webauth(self, ptype = "tool_lab",keyname = "authkey"):
        webcfg = self.config.get_webauth(ptype)
        if keyname in webcfg.keys():
            return webcfg[keyname]
        return ""
        

    def get_mongo_client(self, mtype="default", ref=False, db_version=0):
        try:
            data = web.input()
            if hasattr(data, "db_version"):
                db_version = data.db_version
        except:
            pass
        try:
            db_version = int(db_version)
        except:
            pass
        return self.config.get_mongo_client(mtype, ref, db_version)

    def get_mongo_dbname(self, mtype="default", ref=False, db_version=0):
        try:
            data = web.input()
            if hasattr(data, "db_version"):
                db_version = data.db_version
        except:
            pass
        try:
            db_version = int(db_version)
        except:
            pass
        return self.config.get_mongo_dbname(mtype, ref, db_version)

    def get_biodb_mongo_client(self, mtype="default", db_version=0):
        try:
            data = web.input()
            if hasattr(data, "db_version"):
                db_version = data.db_version
        except:
            pass
        try:
            db_version = int(db_version)
        except:
            pass
        return self.config.get_mongo_client(mtype, True, db_version)

    def get_project_region_bucket(self, project_type="default"):
        return self.config.get_project_region_bucket(project_type)

    def get_bucket_from_path(self, path):
        return self.config.get_bucket_from_path(path)

    def get_rgw_conn(self, region, bucket, new=False):
        return self.config.get_rgw_conn(region, bucket, new)

    def get_work_dir(self):
        return self.config.get_work_dir()


# def get_use_api_clients():
#     return Config().rcf.options("API")
#
#
# def get_api_type(client):
#     if Config().rcf.has_option("API", client):
#         return Config().rcf.get("API", client)
#     else:
#         return None


def get_mongo_client(mtype=None, ref=False):
    return Config().get_mongo_client(mtype, ref)



