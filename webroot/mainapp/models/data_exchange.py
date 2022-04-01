# -*- coding: utf-8 -*-
# __author__ = 'xuting'
import web
import datetime
import json
import re
from mainapp.config.db import Config
# from biocluster.config import Config


class Identity(object):
    """
    验证发送过来的验证码
    """
    def __init__(self, test=False):
        if test:
            self.identity_db = Config().get_identity_db(test=True)
            print "T_IDENTITY_DB", self.identity_db
        else:
            self.identity_db = Config().get_identity_db()
            print "IDENTITY_DB", self.identity_db
        self.record_db = Config().get_record_db()

    def get_task_id(self, code):
        """
        根据验证码， 获取taskid
        """
        where_dict = dict(code=code)
        results = self.identity_db.select("sg_download_code", where=web.db.sqlwhere(where_dict))
        if results:
            for r in results:
                create_time = r["create_time"]
                if (create_time + datetime.timedelta(hours=12)) > datetime.datetime.now():
                    info = {"success": True, "task_id": r["related_task_id"], "info": ""}
                    return info
                else:
                    continue
            info = {"success": False, "task_id": "", "info": "验证码: {} 已经超时， 请重新获取验证码".format(code)}
            return info
        else:
            info = {"success": False, "task_id": "", "info": "验证码: {} 错误！".format(code)}
            return info

    def add_download_record(self, data):
        return self.record_db.insert("download_info", **data)

    def get_target_path(self, code):
        """
        根据验证码，获取要上传的文件的目标路径
        """
        where_dict = dict(code=code)
        results = self.identity_db.select("sg_upload_code", where=web.db.sqlwhere(where_dict))
        if results:
            for r in results:
                create_time = r["create_time"]
                if (create_time + datetime.timedelta(hours=12)) > datetime.datetime.now():
                    info = {"success": True, "rel_path": r["rel_dir"], "info": ""}
                    return info
                else:
                    continue
            info = {"success": False, "rel_path": "", "info": "验证码: {} 已经超时， 请重新获取验证码".format(code)}
            return info
        else:
            info = {"success": False, "rel_path": "", "info": "验证码: {} 错误！".format(code)}
            return info

    def add_upload_record(self, data):
        return self.record_db.insert("upload_info", **data)


class Download(object):
    """
    检索mysql的biocluster库，获取相关的信息， 为分析人员的现在任务文件做准备
    """
    def __init__(self):
        self.table = "workflow"
        self.db = Config().get_db()
        self.client = Config().get_mongo_client()
        self.mongodb = self.client[Config().rcf.get("MONGO", "mongodb")]

    def get_path_by_workflow_id(self, wid):
        """
        更具输入的任务id， 生成相对于workspace的相对路径
        """
        where_dict = dict(workflow_id=wid)
        result = self.db.select(self.table, where=web.db.sqlwhere(where_dict))
        if result:
            r = result[0]
            date = r["run_time"]
            dateStr = "{}{:0>2d}{:0>2d}".format(date.year, date.month, date.day)
            info = json.loads(r["json"])
            name = re.split("\.", info["name"]).pop(-1).split("_")
            name = "".join([i.capitalize() for i in name])
            path = "{}/{}_{}".format(dateStr, name, wid)
            print "相对路径生成完毕，为： {}".format(path)
            return path
        else:
            str_ = "在workflow表里未找到workflow_id为: {} 的记录".format(wid)
            print str_
            return "empty"

    def get_report_by_workflow_id(self, wid):
        collection = self.mongodb['sg_task']
        result = collection.find_one({"task_id": wid})
        if result:
            path = "rerewrweset/files/{}/{}/{}/report_results".format(result['member_id'], result['project_sn'], wid)
            return path
        else:
            return "empty"


if __name__ == "__main__":
    d = Identity()
    d.get_task_id("ASDFGHJKL")
