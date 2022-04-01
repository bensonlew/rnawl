# -*- coding: utf-8 -*-
# __author__ = 'zengjing'

"""
拉取有参RNA demo数据
"""
from biocluster.workflow import Workflow
import pymongo
from biocluster.config import Config
from biocluster.wpm.client import *
import datetime


class RefrnaCopyDemoWorkflow(Workflow):
    """
    last modified by shijin on 20170815
    增加添加备份功能
    """
    def __init__(self, wsheet_object):
        self._sheet = wsheet_object
        super(RefrnaCopyDemoWorkflow, self).__init__(wsheet_object)
        options = [
            {"name": "task_id", "type": 'string', "default": ''},
            {"name": "target_task_id", "type": 'string', "default": ''},
            {"name": "target_project_sn", "type": 'string', "default": ''},
            {"name": "target_member_id", "type": 'string', "default": ''}
        ]
        self.add_option(options)
        self.set_options(self._sheet.options())
        self._project_type = "ref_rna"
        self.client = Config().get_mongo_client(mtype=self._project_type)
        self.db = self.client[Config().get_mongo_dbname(self._project_type)]

    def check(self):
        pass

    def run(self):
        self.start_listener()
        self.fire("start")
        #time = datetime.datetime.now().strftime("%Y-%m-%d %H:%M:%S")
        #old_task_id = self.old_task_id(time)
        #self.update_task_id(old_task_id, self.option("target_task_id"), time)
        #self.logger.info("替换旧task_id: {}完毕".format(old_task_id))
        worker = worker_client()
        id = datetime.datetime.now().strftime("%Y%m%d_%H%M%S%f")[:-3]
        self.logger.info("开始备份新demo，新demo的id为{}".format(self.option("task_id") + '_' + id))
        json_obj = {
            "type": "workflow",
            "id": id,
            "name": "copy_demo.demo_backup",
            "IMPORT_REPORT_DATA": True,
            "IMPORT_REPORT_AFTER_END": False,
            "options": {
                "task_id": self.option("task_id"),
                "target_task_id": self.option("task_id") + '_' + id,
                "target_project_sn": "refrna_demo",
                "target_member_id": "refrna_demo"
            }
        }
        info = worker.add_task(json_obj)
        self.logger.info(info)
        self.end()

    def old_task_id(self, time):
        # db = Config().mongo_client[Config().MONGODB + "_ref_rna"]
        # col = db["sg_task"]
        col = self.db["sg_task"]
        result = col.find_one({"task_id": {"$regex": self.option("task_id") + "_.*_.*"}, "demo_status": "end"})
        if result:
            old_task_id = result["task_id"]
            col.update_one({"_id": result["_id"]}, {"$set": {"task_id": self.option("target_task_id")}})
            col.update_one({"_id": result["_id"]}, {"$set": {"member_id": self.option("target_member_id")}})
            col.update_one({"_id": result["_id"]}, {"$set": {"project_sn": self.option("target_project_sn")}})
            col.update_one({"_id": result["_id"]}, {"$set": {"is_demo": 2}})
            col.update_one({"_id": result["_id"]}, {"$set": {"created_ts": time}})
            return old_task_id
        else:
            raise Exception("没有备份好的demo数据")

    def update_task_id(self, old_task_id, new_task_id, time):
        # db = Config().mongo_client[Config().MONGODB + "_ref_rna"]
        col_list = []
        # for col_name in db.collection_names():
        for col_name in self.db.collection_names():
            col = self.db[col_name]
            result = col.find_one()
            try:
                if "task_id" in result:
                    col_list.append(col_name)
            except:
                continue
        for col_name in col_list:
            col = self.db[col_name]
            results = col.find({"task_id": old_task_id})
            for result in results:
                col.update_one({"_id": result["_id"]}, {"$set": {"task_id": new_task_id}})
                col.update_one({"_id": result["_id"]}, {"$set": {"created_ts": time}})

    def end(self):
        super(RefrnaCopyDemoWorkflow, self).end()
