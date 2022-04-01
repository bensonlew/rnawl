# !usr/bin/python
# -*- coding: utf-8 -*-
# __author__ = 'WangWenjie'

import os
import re
import time
import shutil
import json
import unittest
import glob
import datetime
import csv
from biocluster.workflow import Workflow
from biocluster.api.file.lib.transfer import MultiFileTransfer
from biocluster.core.exceptions import OptionError
from biocluster.config import Config
from bson.objectid import ObjectId

class LmplotWorkflow(Workflow):
    """
    标准曲线流程
    """
    def __init__(self, wsheet_object):
        self._sheet = wsheet_object
        super(LmplotWorkflow, self).__init__(wsheet_object)
        options = [
            {"name":"main_id", "type":"string"},
            {"name":"list_file","type":"infile","format":"denovo_rna_v2.common_dir"},
            {"name":"update_info", "type": "string"},
        ]
        self.add_option(options)
        self.set_options(self._sheet.options())
        # self.lmplot = self.add_tool("basic_molecular.lmplot")

    def check_option(self):
        """
        参数检查
        """
        if not self.option("list_file"):
            raise OptionError("没有基因文件输入")  

    def run_lmplot(self):
        """
        获取line.txt文件和gene_name,进行标准曲线计算
        """       
        list_file = os.path.join(self.option("list_file").prop["path"])
        self.lmplot_tools = []
        filelist = os.listdir(list_file)
        for f in filelist:
            f = os.path.join(list_file, f)
            dirname = os.path.dirname(f)
            baseName = os.path.basename(f)
            if dirname.endswith(os.sep):
                options = {
                    "input_file": dirname+baseName,
                    "gene_name": baseName.split('.')[0]
                }
            else:
                options = {
                    "input_file": dirname+os.sep+baseName,
                    "gene_name": baseName.split('.')[0]
                }
            lmplot_tool = self.add_tool("basic_molecular.lmplot")
            lmplot_tool.set_options(options)
            self.lmplot_tools.append(lmplot_tool)
        if len(self.lmplot_tools) == 1:
            self.lmplot_tools[0].on('end',self.set_output)
        else:
            self.on_rely(self.lmplot_tools,self.set_output)
        for tool in self.lmplot_tools:
            tool.run()

    def set_output(self,event):
        self.logger.info("设置结果目录")
        os.mkdir(self.output_dir + "/lmplot")
        for lmplot in self.lmplot_tools:
            files =os.listdir(lmplot.output_dir)
            for file in files:
                print file
                if os.path.exists(self.output_dir + '/' + file):
                    os.remove(self.output_dir + '/' + file)
                os.link(lmplot.output_dir + '/' + file, self.output_dir + '/lmplot/' + file)
        os.mkdir(self.output_dir + "/png_dir")
        os.mkdir(self.output_dir + "/txt_dir")
        os.system('cd {} && cp *.png {}'.format(self.output_dir + "/lmplot", self.output_dir + "/png_dir"))
        os.system('cd {} && cp *.txt {}'.format(self.output_dir + "/lmplot", self.output_dir + "/txt_dir"))
        self.set_db()

    def set_db(self):
        self.logger.info("开始导表")
        api_lmplot = self.api.api("basic_molecular.lmplot")
        s3_upload_path = '{}{}'.format(self._sheet.output, 'png_dir/')
        api_lmplot.update_lmplot_record(self.option('main_id'), self.output_dir + '/txt_dir', s3_upload_path)
        self.logger.info("导表结束")
        self.end()

    def run(self):
        self.run_lmplot()
        super(LmplotWorkflow, self).run()

    def end(self):
        result_dir = self.add_upload_dir(self.output_dir)
        result_dir.add_relpath_rules([
            [".", "", "结果输出目录"],
        ])
        result_dir.add_regexp_rules([
            ["", "", ""]
        ])
        super(LmplotWorkflow, self).end()
