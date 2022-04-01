# -*- coding: utf-8 -*-
# __author__ = 'zhaozhigang'

import os
from biocluster.workflow import Workflow
import datetime
import unittest
import types
import shutil
from bson.objectid import ObjectId


class TrinityWorkflow(Workflow):
    """
    基于reads的humann2分析

    """

    def __init__(self, wsheet_object):
        self._sheet = wsheet_object
        super(TrinityWorkflow, self).__init__(wsheet_object)
        options = [
            {'name': 'in_fastq', 'type': 'infile', 'format': 'sequence.fastq_dir'},  # 输入的fq序列文件夹
            {'name': 'min_contig_length', 'type': 'int', 'default': 200},  # trnity报告出的最短的Contig长度，默认为200
            {'name': 'kmer_size', "type": 'int', 'default': 25},  # Trinity软件使用的kmer的长度，默认为25
            {'name': 'min_kmer_cov', 'type': 'int', 'default': 2},  # 最小kmer的计数
            {"name": "main_id", "type": "string"},
            {"name": "update_info", "type": "string"},
        ]
        self.add_option(options)
        self.revise_infiles()
        self.set_options(self._sheet.options())

    def run(self):
        self.run_main()
        super(TrinityWorkflow, self).run()

    def run_main(self):
        self.main = self.add_module("toolapps.trinity")
        options = {
            "in_fastq": self.option('in_fastq'),
            "min_contig_length": self.option('min_contig_length'),
            "kmer_size": self.option('kmer_size'),
            "min_kmer_cov": self.option('min_kmer_cov'),
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
        api_main.add_main_table("trinity", main_id = self.option('main_id'))
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
            [".", "", "Trinity组装拼接结果文件目录",0],
        ])
        super(TrinityWorkflow, self).end()