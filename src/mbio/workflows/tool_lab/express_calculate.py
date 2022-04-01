# -*- coding: utf-8 -*-
# __author__ = 'zhaozhigang'

import os
from biocluster.workflow import Workflow
import datetime
import unittest
import types
import shutil
from bson.objectid import ObjectId


class ExpressCalculateWorkflow(Workflow):
    """
    基于reads的humann2分析

    """

    def __init__(self, wsheet_object):
        self._sheet = wsheet_object
        super(ExpressCalculateWorkflow, self).__init__(wsheet_object)
        options = [
            {"name": "fq_type", "type": "string"},  # PE OR SE
            {"name": "ref_fa", "type": "infile", "format": "sequence.fasta"},  # trinity.fasta文件
            {"name": "ref_dir", "type": "infile", "format": "toolapps.fasta_dir"},  # trinity.fasta文件夹
            {"name": "fastq_dir", "type": "infile", "format": "sequence.fastq_dir"}, # PE数据或者SE数据，包含所有的左端和右端的fq文件或者single端文件
            {"name": "exp_way", "type": "string", "default": "fpkm"},  # 计算表达量的指标或者方法
            {"name": "software", "type": "string", 'default': "rsem"},  ## 计算使用的软件
            {"name": "main_id", "type": "string"},
            {"name": "update_info", "type": "string"},
        ]
        self.add_option(options)
        self.revise_infiles()
        self.set_options(self._sheet.options())

    def run(self):
        self.run_main()
        super(ExpressCalculateWorkflow, self).run()

    def run_main(self):
        self.main = self.add_module("toolapps.express_calculate")
        options = {
            "fq_type": self.option('fq_type'),
            "ref_fa": self.option('ref_fa'),
            "ref_dir": self.option('ref_dir'),
            "fastq_dir": self.option('fastq_dir'),
            "exp_way": self.option('exp_way'),
            "software": self.option('software'),
        }
        self.main.set_options(options)
        self.main.on('end', self.set_output)
        self.main.run()

    def set_output(self):
        self.linkdir(self.main.output_dir)
        self.set_db()

    def set_db(self):
        self.logger.info("开始导表")
        api_main = self.api.api("tool_lab.common")
        api_main.add_main_table("express_calculate", main_id = self.option('main_id'))
        self.logger.info("导表完成")
        self.end()

    def linkdir(self, dirpath):
        allfiles = os.listdir(dirpath)
        newdir = self.output_dir
        if not os.path.exists(newdir):
            os.mkdir(newdir)
        oldfiles = [os.path.join(dirpath, i) for i in allfiles]
        newfiles = [os.path.join(newdir, i) for i in allfiles]
        for newfile in newfiles:
            if os.path.exists(newfile):
                if os.path.isfile(newfile):
                    os.remove(newfile)
                else:
                    os.system('rm -r %s' % newfile)
                    # self.logger.info('rm -r %s' % newfile)
        for i in range(len(allfiles)):
            if os.path.isfile(oldfiles[i]):
                os.link(oldfiles[i], newfiles[i])
            elif os.path.isdir(oldfiles[i]):
                # self.logger.info('cp -r %s %s' % (oldfiles[i], newdir))
                os.system('cp -r %s %s' % (oldfiles[i], newdir))

    def end(self):
        result_dir = self.add_upload_dir(self.output_dir)
        self.logger.info("开始上传结果文件")
        result_dir.add_relpath_rules([
            [".", "", "表达量分析模块结果输出目录",0],
        ])
        super(ExpressCalculateWorkflow, self).end()