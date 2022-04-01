# -*- coding: utf-8 -*-
# __author__ = 'zengjing'
# __last_modified__ = 20180227
from biocluster.workflow import Workflow
import pymongo
from biocluster.config import Config
from biocluster.wpm.client import *
import datetime


class DemoInitWorkflow(Workflow):
    """
    初始化设置demo时进行demo备份
    type：ref_rna时，对ref_rna项目做初始化，进行备份或者删除备份；
          ref_rna_demo时，删除拉取的ref_rna任务（前端做计划任务，传递拉取已有30天的demo参数）
    """
    def __init__(self, wsheet_object):
        self._sheet = wsheet_object
        super(DemoInitWorkflow, self).__init__(wsheet_object)
        options = [
            {"name": "task_id", "type": "string"},  # 要设置为demo或取消的demo的task_id
            {"name": "type", "type": "string", "default": "ref_rna"},  # demo的类型
            {"name": "setup_type", "type": "string", "default": "setup"},  # 对demo进行的操作，设置为demo，取消删除demo
            {"name": "demo_number", "type": "int", "default": 30}  # demo备份的数量
        ]
        self.add_option(options)
        self.set_options(self._sheet.options())

    def check_options(self):
        if self.option("task_id") == "" or self.option("task_id") == " ":
            raise OptionError("task_id不能为空")

    def run(self):
        self.start_listener()
        self.fire("start")
        if self.option("type") == "ref_rna":
            if self.option("setup_type") == "setup":
                self._project_type = "ref_rna"
                self.client = Config().get_mongo_client(mtype=self._project_type)
                self.db = self.client[Config().get_mongo_dbname(self._project_type)]
                from mbio.packages.rna.refrna_copy_demo import RefrnaCopyMongo
                target_project_sn = "refrna_demo"
                target_member_id = "refrna_demo"
                self.logger.info("开始备份新的demo")
                for i in range(self.option("demo_number")):
                    target_task_id = self.option("task_id") + "_" + datetime.datetime.now().strftime("%Y%m%d_%H%M%S%f")[:-3]
                    copy_task = RefrnaCopyMongo(self.option("task_id"), target_task_id, target_project_sn, target_member_id)
                    copy_task.run()
                    # db = Config().mongo_client[Config().MONGODB + "_ref_rna"]
                    col = self.db["sg_task"]
                    result = col.find_one({"task_id": target_task_id, "project_sn": target_project_sn})
                    col.update_one({"_id": result["_id"]}, {"$set": {"demo_status": "end"}})
                self.logger.info("备份{}份新的demo成功".format(self.option("demo_number")))
            if self.option("setup_type") in ["cancel", "delete"]:
                self.logger.info("开始删除备份的demo")
                from mbio.packages.rna.refrna_copy_delete import RefrnaCopyDelete
                RefrnaCopyDelete().find_task_id(task_id=self.option("task_id"))
        if self.option("type") == "ref_rna_demo" and self.option("setup_type") == "delete":
            self.logger.info("开始删除拉取的demo")
            from mbio.packages.rna.refrna_copy_delete import RefrnaCopyDelete
            RefrnaCopyDelete().remove(self.option("task_id"))
        if self.option("type") == "metagenomic_demo" and self.option("setup_type") == "delete":
            self.logger.info("开始删除拉取的demo")
            from mbio.packages.metagenomic.delete_demo import DeleteDemo
            DeleteDemo().remove(self.option("task_id"))  # add metagenomic_demo by GHD @20180227
        self.end()

    def end(self):
        super(DemoInitWorkflow, self).end()
