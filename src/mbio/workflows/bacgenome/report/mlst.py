# -*- coding: utf-8 -*-
# __author__ = 'hao.gao'

import os
import re
import types
import shutil
from biocluster.workflow import Workflow
from bson import ObjectId
import gevent
from mbio.packages.bacgenome.common import link_dir,link_file


class MlstWorkflow(Workflow):
    """
    细菌基因组mlst分析
    """
    def __init__(self, wsheet_object):
        self._sheet = wsheet_object
        self.rpc = False
        super(MlstWorkflow, self).__init__(wsheet_object)
        options = [
            {"name": "specimen_id", "type": "string"},  # 样品名称
            {"name": "species", "type": "string", },  # 物种名称
            {"name": "task_id", "type": "string"},
            {"name": "update_info", "type": "string"},
            {"name": "main_id", "type": "string"},
        ]
        self.add_option(options)
        self.set_options(self._sheet.options())
        self.common_api = self.api.api('bacgenome.common_api')
        self.modules = []

    def download_file(self):
        """
        download file from s3
        :return:
        """
        self.logger.info("开始下载组装序列文件！")
        assemble_dir = os.path.join(self.work_dir, "assemble_dir")
        if os.path.exists(assemble_dir):
            shutil.rmtree(assemble_dir)
        self.samples = str(self.option("specimen_id")).split(",")
        (self.assemble_dir, self.analysis_type)= self.common_api.down_seq_files(assemble_dir, self.option("task_id"), self.samples)
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
        for sample in self.samples:
            sample_path = os.path.join(assemble_dir2, sample + "_scaf.fna")
            assemble_dir = os.path.join(self.assemble_dir, sample)
            dir_list = os.listdir(assemble_dir)
            for file2 in dir_list:
                os.system("cat {} >> {}".format(os.path.join(assemble_dir, file2), sample_path))

    def run_mlst(self):
        """
        MLST预测
        :return:
        """
        self.logger.info("开始用mlst软件进行预测！")
        self.samples = str(self.option("specimen_id")).split(",")
        self.download_file()
        for sample in self.samples:
            self.mlst = self.add_tool('bacgenome.mlst')
            if self.analysis_type in ['uncomplete']:
                fa = os.path.join(self.work_dir, "assemble_dir", sample +"_scaf.fna")
            else:
                fa = os.path.join(self.work_dir, "assemble_dir2", sample + "_scaf.fna")
            options = {
                "fasta": fa,
                "species": self.option("species"),
                "sample_name": sample
                }
            self.mlst.set_options(options)
            self.modules.append(self.mlst)
        if len(self.modules) > 1:
            self.on_rely(self.modules, self.set_output)
        else:
            self.modules[0].on('end', self.set_output)
        for module in self.modules:
            module.run()
            gevent.sleep(0)

    def run(self):
        self.IMPORT_REPORT_DATA = True
        self.IMPORT_REPORT_AFTER_END = False
        self.run_mlst()
        super(MlstWorkflow, self).run()

    def set_output(self):
        """
        设置结果目录文件
        :return:
        """
        self.logger.info("开始设置结果文件目录！")
        for module in self.modules:
            for i in os.listdir(module.output_dir):
                link_file(module.output_dir + "/" + i, self.output_dir + "/" + i)
        self.logger.info("设置结果文件目录完成!")
        self.set_db()

    def set_db(self):
        """
        导入MongoDB数据
        :return:
        """
        self.logger.info("start mongo>>>>>>>>>>>>")
        api_mlst = self.api.api('bacgenome.mlst')
        main_id = self.option("main_id")
        for sample in self.samples:
            if os.path.exists(self.output_dir + "/" +sample + ".mlst.ST.xls"):
                api_mlst.add_mlst_stat(main_id, self.output_dir + "/" + sample + ".mlst.ST.xls")
                api_mlst.add_mlst_detail(main_id, self.output_dir + "/" + sample + ".mlst.detail.xls", sample)
        self.logger.info("end MongoDB<<<<<<<<<<<<<")
        self.end()

    def end(self):
        repaths = [
            [".", "", "MLST分析结果目录",0],
        ]
        regexps = [
            [r'.*.mlst.ST.xls', 'xls', 'MLST分析统计表',0],
            [r'.*.mlst.detail.xls', 'xls', 'MLST分析详情表',0],
            [r'.*.HitMLST.fasta', 'xls', 'MLST分析的基因序列',0],
        ]
        sdir = self.add_upload_dir(self.output_dir)
        sdir.add_relpath_rules(repaths)
        sdir.add_regexp_rules(regexps)
        super(MlstWorkflow, self).end()