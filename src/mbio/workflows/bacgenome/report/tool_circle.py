# -*- coding: utf-8 -*-
# __author__ = 'hao.gao'

import os
import re
from bson import SON
import shutil
from biocluster.workflow import Workflow
from bson import ObjectId
import datetime
import gevent
from mbio.packages.bacgenome.common import link_dir,link_file
from biocluster.config import Config
import HTMLParser

class ToolCircleWorkflow(Workflow):
    """
    和弦图细菌基因组COG/GO/KEGG注释
    """
    def __init__(self, wsheet_object):
        self._sheet = wsheet_object
        super(ToolCircleWorkflow, self).__init__(wsheet_object)
        options = [
            {"name": "task_id", "type": "string"},
            {"name": "update_info", "type": "string"},
            {"name": "main_id", "type": "string"},
            {"name": "specimens", "type": "string"},
            {"name": "data_type", "type": "string"},
            {"name": "relate_name", "type": "string"},
            {"name": "params", "type": "string"},
            {"name": "tool_type", "type": "string"},
        ]
        self.add_option(options)
        self.set_options(self._sheet.options())
        self.common_api = self.api.api('bacgenome.common_api')
        self.modules = []
        self.file_path = self._sheet.output
        self.samples = {}

    def download_file(self):
        """
        download file from s3
        :return:
        """
        self.logger.info("开始下载组装序列文件！")
        anno_dir = os.path.join(self.work_dir, "anno_dir")
        if os.path.exists(anno_dir):
            shutil.rmtree(anno_dir)
        samples = eval(HTMLParser.HTMLParser().unescape(self.option("specimens")))
        for i in samples:
            self.samples[i['id']] =i['name']
        self.anno_dir = self.common_api.get_anno_file(anno_dir, self.option("task_id"),self.samples , self.option("data_type"))


    def run_circle(self):
        """
        :return:
        """
        self.logger.info("开始用预测！")
        self.download_file()
        self.tool_circle = self.add_tool('bacgenome.tool_circle')
        options = {
            "file_dir": self.anno_dir,
            "data_type": self.option("data_type")
        }
        self.tool_circle.set_options(options)
        self.tool_circle.on("end", self.set_output)
        self.tool_circle.run()

    def run(self):
        self.IMPORT_REPORT_DATA = True
        self.IMPORT_REPORT_DATA_AFTER_END = False
        self.run_circle()
        super(ToolCircleWorkflow, self).run()

    def set_output(self):
        """
        设置结果目录文件
        :return:
        """
        self.logger.info("开始设置结果文件目录！")
        if os.path.exists(self.output_dir+"/" + "all.table.xls"):
            os.remove(self.output_dir+"/" + "all.table.xls")
        link_file(self.tool_circle.output_dir +"/" + "all.table.xls", self.output_dir+"/" + "all.table.xls")
        self.remot1 = self.file_path + '/' + "all.table.xls"
        self.logger.info("设置结果文件目录完成!")
        self.set_db()

    def set_db(self):
        """
        导入MongoDB数据
        :return:
        """
        self.logger.info("start mongo>>>>>>>>>>>>")
        dict= {'path': self.remot1}
        self.common_api.update_data(self.option("main_id"), "sg_tool_lab_circle", dict)
        self.logger.info("end MongoDB<<<<<<<<<<<<<")
        self.end()

    def end(self):
        self.add_upload_dir(self.output_dir)
        try:
            main_id = ObjectId(self.option('main_id'))
        except:
            pass
        main_info = dict(
            task_id=self.option('task_id'),
            project_type='bacgenome',
            params=self.option('params'),
            status="end",
            main_id=main_id,
            tool_type=self.option('tool_type'),
            relate_id=main_id,
            relate_name=self.option('relate_name'),
            file_path=self.remot1
        )
        collection_name = 'tool_thurl'
        project_type = 'tool_lab'
        db = Config().get_mongo_client(mtype=project_type)[Config().get_mongo_dbname(project_type)]
        main_id_ = db[collection_name].insert_one(SON(main_info)).inserted_id
        conn = db[collection_name]
        task_id = conn.find_one({'_id': main_id_})['task_id']
        conn.update({'_id': main_id_, "task_id": task_id}, {"$set": {'main_id': main_id_}}, upsert=True)
        super(ToolCircleWorkflow, self).end()