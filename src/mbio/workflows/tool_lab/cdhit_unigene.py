# -*- coding: utf-8 -*-
# __author__ = 'zhaozhigang'

import os
from biocluster.workflow import Workflow
import datetime
import unittest
import types
import shutil
from bson.objectid import ObjectId


class CdhitUnigeneWorkflow(Workflow):
    """
    基于reads的humann2分析

    """

    def __init__(self, wsheet_object):
        self._sheet = wsheet_object
        super(CdhitUnigeneWorkflow, self).__init__(wsheet_object)
        options = [
            {"name": "gene_fa_dir", "type": "infile", "format": "toolapps.fasta_dir"},  # 基因预测核酸文件
            {"name": "gene_faa_dir", "type": "infile", "format": "toolapps.fasta_dir"},  # 基因预测蛋白文件
            {"name": "cdhit_identity", "type": "float", "default": 0.9},  # 给出cdhit的参数identity
            {"name": "cdhit_coverage", "type": "float", "default": 0.9},  # 给出cdhit的参数coverage
            {"name": "ana_type", "type": "string", "default": "nucl"},  # 输入分析类型，对蛋白还是对核酸聚类
            {"name": "main_id", "type": "string"},
            {"name": "update_info", "type": "string"},
        ]
        self.add_option(options)
        self.revise_infiles()
        self.set_options(self._sheet.options())

    def run(self):
        self.run_main()
        super(CdhitUnigeneWorkflow, self).run()

    def run_main(self):
        self.main = self.add_module("toolapps.cdhit_unigene")
        options = {
            "gene_fa_dir": self.option('gene_fa_dir'),
            "gene_faa_dir": self.option('gene_faa_dir'),
            "cdhit_identity": self.option('cdhit_identity'),
            "cdhit_coverage": self.option('cdhit_coverage'),
            "ana_type": self.option('ana_type'),
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
        api_main.add_main_table("cdhit_unigene", main_id = self.option('main_id'))
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
        self.logger.info("Now start upload!")
        result_dir.add_relpath_rules([
            [".", "", "非冗余基因集输出目录",0],
            ["./geneCatalog_stat.xls", "xls", "非冗余基因集统计结果",0],
            ["./gene.uniGeneset.fa", "fa", "非冗余基因集核酸序列",0],
            ["./gene.uniGeneset.faa", "faa", "非冗余基因集蛋白序列",0],
        ])
        super(CdhitUnigeneWorkflow, self).end()