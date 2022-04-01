# -*- coding: utf-8 -*-
# __author__ = 'qingchen.zhang'


import os
import re
import types
from biocluster.workflow import Workflow
from bson import ObjectId
import datetime
import gevent
import shutil
from mbio.packages.bacgenome.common import link_dir


class ResfinderWorkflow(Workflow):
    """
    细菌基因组resfinder耐药基因预测
    """
    def __init__(self, wsheet_object):
        self._sheet = wsheet_object
        self.rpc = False
        super(ResfinderWorkflow, self).__init__(wsheet_object)
        options = [
            {"name": "specimen_id", "type": "string"},  # 样品名称
            {"name": "task_id", "type": "string"},
            {"name": "update_info", "type": "string"},
            {"name": "main_id", "type": "string", 'default': ""},
            {"name": "coverage", "type": "float", "default": 60.0},# 最小的coverage
            {"name": "identity", "type": "float", "default": 80.0},## 最小的identity
            {"name": "species_name", "type": "string", "default": ""}, ## 点突变物种名称
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
        self.logger.info("开始下载基因序列文件！")
        gene_dir = os.path.join(self.work_dir, "gene_dir")
        if os.path.exists(gene_dir):
            shutil.rmtree(gene_dir)
        self.gene_dir = self.common_api.down_gene_files(gene_dir, self.option("task_id"), self.samples)

    def run_resfinder(self):
        """
        resfinder预测
        :return:
        """
        self.logger.info("开始用resfinder软件进行预测！")
        self.samples = str(self.option("specimen_id")).split(",")
        self.logger.info(self.samples)
        self.download_file()
        if self.option("species_name") in ["", "All", 'all', "ALL"]:
            species_name = ""
        else:
            species_name = self.option("species_name")
        for sample in self.samples:
            self.resfinder_predict = self.add_tool('bacgenome.resfinder')
            options = {
                "gene_fa": os.path.join(self.gene_dir, sample, sample + ".fnn"),
                "gene_gff": os.path.join(self.gene_dir, sample, sample + ".gff"),
                "sample_name": sample,
                "min_cov": self.option("coverage")/100,
                "min_iden": self.option("identity")/100,
                "resfinder_database": "true",
                "disinfinder_database": "true",
                "species_name": species_name
            }
            self.resfinder_predict.set_options(options)
            self.modules.append(self.resfinder_predict)
        if len(self.modules) > 1:
            self.on_rely(self.modules, self.set_output)
        else:
            self.modules[0].on('end', self.set_output)
        for module in self.modules:
            module.run()
            gevent.sleep(0)

    def run(self):
        self.IMPORT_REPORT_DATA = True
        self.IMPORT_REPORT_DATA_AFTER_END = False
        self.run_resfinder()
        super(ResfinderWorkflow, self).run()

    def set_output(self):
        """
        设置结果目录文件
        :return:
        """
        self.logger.info("开始设置结果文件目录！")
        for module in self.modules:
            for sample_dir in os.listdir(module.output_dir):
                sample_dir_path = module.output_dir + "/" +sample_dir
                output_sample = self.output_dir + "/" +sample_dir
                if os.path.exists(output_sample):
                    shutil.rmtree(output_sample)
                if os.path.exists(sample_dir_path):
                    if len(os.listdir(sample_dir_path)) != 0:
                        link_dir(sample_dir_path, output_sample)
        self.logger.info("设置结果文件目录完成!")
        self.set_db()


    def set_db(self):
        self.logger.info("start mongo>>>>>>>>>>>>")
        api_resfinder = self.api.api('bacgenome.resfinder')
        if self.option("main_id") not in [""]:
            main_id = self.option("main_id")
        else:
            params = {
                "group_detail": {"All": self.samples},
                "identity": float(80),
                "coverage": float(60),
                "species_name": "All"
            }
            main_id = api_resfinder.add_resfinder(params=params)
        for sample in self.samples:
            self.logger.info(sample)
            sample_dir = os.path.join(self.output_dir, sample)
            api_resfinder.add_resfinder_dir(sample_dir, sample, main_id=main_id)
        self.end()

    def end(self):
        repaths = [
            [".", "", "耐药基因Resfinder预测结果目录",0,""],
        ]
        regexps = [
            [r'.*_stat.xls', 'xls', '耐药基因Resfinder预测结果统计表',0,""],
            [r'.*_resfinder.class.xls', 'xls', 'Resfinder分类统计表',0,""],
            [r'.*_disinfinder.class.xls', 'xls', 'Desinfinder分类统计表',0,""],
            [r'.*_resfinder.detail.xls', 'xls', 'Resfinder分类详情表',0,""],
            [r'.*_disinfinder.detail.xls', 'xls', 'Desinfinder分类详情表',0,""]
        ]
        sdir = self.add_upload_dir(self.output_dir)
        sdir.add_relpath_rules(repaths)
        sdir.add_regexp_rules(regexps)
        super(ResfinderWorkflow, self).end()
