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


class ToolNglocWorkflow(Workflow):
    """
    细菌基因组打通的原核亚细胞定位ngloc
    """
    def __init__(self, wsheet_object):
        self._sheet = wsheet_object
        super(ToolNglocWorkflow, self).__init__(wsheet_object)
        options = [
            {"name": "specimen_id", "type": "string"},  # 样品名称
            {"name": "task_id", "type": "string"},
            {"name": "update_info", "type": "string"},
            {"name": "main_id", "type": "string"},
            {"name": "specimen_id_new", "type": "string"},
            {"name": "relate_name", "type": "string"},
            {"name": "params", "type": "string"},
            {"name": "tool_type", "type": "string"},
        ]
        self.add_option(options)
        self.set_options(self._sheet.options())
        self.common_api = self.api.api('bacgenome.common_api')
        self.modules = []
        self.file_path = self._sheet.output

    def download_file(self):
        """
        download file from s3
        :return:
        """
        self.logger.info("开始下载基因序列文件！")
        (self.faa_file, self.pfam_file)= self.common_api.get_gene_ngloc(self.work_dir, self.option("task_id"), self.option("specimen_id"))

    def run_ngloc(self):
        """
        原核亚细胞定位ngloc文件处理
        :return:
        """
        self.logger.info("开始用预测！")
        self.sample = self.option("specimen_id_new")
        self.download_file()
        self.tool_ngloc = self.add_tool('bacgenome.tool_ngloc')
        options = {
            "gene_faa": self.faa_file,
            "pfam": self.pfam_file,
            'sample': self.sample,
            }
        self.tool_ngloc.set_options(options)
        self.tool_ngloc.on("end", self.set_output)
        self.tool_ngloc.run()

    def run(self):
        self.IMPORT_REPORT_DATA = True
        self.IMPORT_REPORT_DATA_AFTER_END = False
        self.run_ngloc()
        super(ToolNglocWorkflow, self).run()

    def set_output(self):
        """
        设置结果目录文件
        :return:
        """
        self.logger.info("开始设置结果文件目录！")
        link_file(self.tool_ngloc.output_dir + '/' + self.sample+".pfam.fasta", self.output_dir+ '/' + self.sample+".pfam.fasta")
        self.remot =self.file_path + '/' + self.sample+".pfam.fasta"
        self.logger.info("设置结果文件目录完成!")
        self.set_db()

    def set_db(self):
        """
        导入MongoDB数据
        :return:
        """
        self.logger.info("start mongo>>>>>>>>>>>>")
        dict= {'path': self.remot}
        self.common_api.update_data(self.option("main_id"), "sg_tool_lab_ngloc", dict)
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
            file_path=self.remot
        )
        collection_name = 'tool_thurl'
        project_type = 'tool_lab'
        db = Config().get_mongo_client(mtype=project_type)[Config().get_mongo_dbname(project_type)]
        main_id_ = db[collection_name].insert_one(SON(main_info)).inserted_id
        conn = db[collection_name]
        task_id = conn.find_one({'_id': main_id_})['task_id']
        conn.update({'_id': main_id_, "task_id": task_id}, {"$set": {'main_id': main_id_}}, upsert=True)
        super(ToolNglocWorkflow, self).end()