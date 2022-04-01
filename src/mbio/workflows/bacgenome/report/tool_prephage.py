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


class ToolPrephageWorkflow(Workflow):
    """
    细菌基因组打通的引物设计
    """
    def __init__(self, wsheet_object):
        self._sheet = wsheet_object
        super(ToolPrephageWorkflow, self).__init__(wsheet_object)
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
        self.logger.info("开始下载组装序列文件！")
        assemble_dir = os.path.join(self.work_dir, "assemble_dir")
        if os.path.exists(assemble_dir):
            shutil.rmtree(assemble_dir)
        (self.assemble_dir, self.analysis_type)= self.common_api.down_seq_files(assemble_dir, self.option("task_id"), [self.option("specimen_id")])
        if self.analysis_type in ['complete']:
            self.cat_reads()

    def cat_reads(self):
        """
        合并reads 完成图合并fasta文件，形成单个样本文件一个
        :return:
        """
        assemble_dir2 = os.path.join(self.work_dir, "assemble_dir2")
        if os.path.exists(assemble_dir2):
            shutil.rmtree(assemble_dir2)
        os.mkdir(assemble_dir2)
        for sample in [self.option("specimen_id")]:
            sample_path = os.path.join(assemble_dir2, sample + "_scaf.fna")
            assemble_dir = os.path.join(self.assemble_dir, sample)
            dir_list = os.listdir(assemble_dir)
            for file2 in dir_list:
                os.system("cat {} >> {}".format(os.path.join(assemble_dir, file2), sample_path))

    def run(self):
        self.IMPORT_REPORT_DATA = True
        self.IMPORT_REPORT_DATA_AFTER_END = False
        self.start_listener()
        self.fire("start")
        self.set_output()
        super(ToolPrephageWorkflow, self).run()

    def set_output(self):
        """
        设置结果目录文件
        :return:
        """
        self.logger.info("开始用预测！")
        self.samples = [self.option("specimen_id")]
        self.download_file()
        self.table = self.common_api.get_table_file(self.work_dir, self.option("task_id"), self.option("specimen_id"),
                                                    "prephage")
        with open(self.table) as f, open(self.work_dir + "/" + self.option('specimen_id_new') + ".prephage.xls",
                                         "w") as g:
            g.write("Location\tElement ID\tStart\tEnd\tLength\n")
            lines = f.readlines()
            for line in lines[1:]:
                lin = line.strip().split("\t")
                g.write("{}\t{}\t{}\t{}\t{}\n".format(lin[1], lin[0], lin[2], lin[3], lin[4]))
        self.logger.info("开始设置结果文件目录！")
        link_file(self.work_dir+"/"+self.option('specimen_id_new') + ".prephage.xls", self.output_dir+"/"+self.option('specimen_id_new') + ".prephage.xls")
        if os.path.exists(self.output_dir + "/" + self.option('specimen_id_new') + ".fasta"):
            os.remove(self.output_dir + "/" + self.option('specimen_id_new') + ".fasta")
        if self.analysis_type in ['uncomplete']:
            os.link(os.path.join(self.assemble_dir, self.option("specimen_id") + "_scaf.fna"), self.output_dir+"/"+self.option('specimen_id_new')+".fasta")
        else:
            os.link(os.path.join(self.work_dir,"assemble_dir2", self.option("specimen_id") + "_scaf.fna"),
                    self.output_dir + "/" + self.option('specimen_id_new') + ".fasta")
        self.remot1 =self.file_path +'/'+ self.option('specimen_id_new')+".fasta"
        self.remot2 = self.file_path + '/' + self.option('specimen_id_new') + ".prephage.xls"
        self.logger.info("设置结果文件目录完成!")
        self.set_db()

    def set_db(self):
        """
        导入MongoDB数据
        :return:
        """
        self.logger.info("start mongo>>>>>>>>>>>>")
        dict= {'path': self.remot1,'path1':self.remot2}
        self.common_api.update_data(self.option("main_id"), "sg_tool_lab_prephage", dict)
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
            path=self.remot1,
            path1=self.remot2,
        )
        collection_name = 'tool_thurl'
        project_type = 'tool_lab'
        db = Config().get_mongo_client(mtype=project_type)[Config().get_mongo_dbname(project_type)]
        main_id_ = db[collection_name].insert_one(SON(main_info)).inserted_id
        conn = db[collection_name]
        task_id = conn.find_one({'_id': main_id_})['task_id']
        conn.update({'_id': main_id_, "task_id": task_id}, {"$set": {'main_id': main_id_}}, upsert=True)
        super(ToolPrephageWorkflow, self).end()