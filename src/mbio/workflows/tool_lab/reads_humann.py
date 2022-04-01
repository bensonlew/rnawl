# -*- coding: utf-8 -*-
# __author__ = 'zhaozhigang'

import os
from biocluster.workflow import Workflow
import datetime
import unittest
import types
import shutil
from bson.objectid import ObjectId


class ReadsHumannWorkflow(Workflow):
    """
    基于reads的humann2分析

    """

    def __init__(self, wsheet_object):
        self._sheet = wsheet_object
        super(ReadsHumannWorkflow, self).__init__(wsheet_object)
        options = [
            {'name': 'in_fastq', 'type': 'infile', 'format': 'sequence.fastq_dir'},  # 输入的fq文件夹
            {'name': 'qc', 'type': 'bool', 'default': False},  # 是否需要质控
            {'name': 'qc_quality', 'type': 'int', 'default': 20},  # 质控质量值标准
            {'name': 'qc_length', 'type': 'int', 'default': 30},  # 质控最短序列长度
            {'name': 'rm_host', 'type': 'bool', 'default': False},  # 是否需要去除宿主
            {'name': 'ref_database', 'type': 'string', 'default': ''},  # 宿主参考序列库中对应的物种名，eg：E.coli ,B.taurus
            {'name': 'second_ref', 'type': 'string', 'default': ''},  # 参考数据库二级下拉框数据
            {'name': 'ref_undefined', "type": 'infile', 'format': 'sequence.fasta_dir'},
            {'name': 'ref_undefined_name', 'type': 'string', 'default': 'undefined'},  # 自定义参考宿主名称，适应页面参数
            {'name': 'clean_fq', 'type': 'infile', 'format': 'sequence.fastq_dir'},  # 输入的fq文件夹
            {"name": "search_mode", "type": "string", "default": "uniref50"},  # uniref50 or uniref90
            {"name": "prescreen_threshold", "type": 'float', "default": 0.01},
            {"name": "identity_threshold", "type": 'int', "default": 50},
            {"name": "subject_coverage", "type": 'int', "default": 50},
            {"name": "query_coverage", "type": 'int', "default": 90},
            {"name": "pathways", "type": "string", "default": "metacyc"},  # metacyc,unipathway
            {"name": "aligin", "type": "string", "default": "diamond"},  # usearch,rapsearch,diamond
            {"name": "main_id", "type": "string"},
            {"name": "update_info", "type": "string"},
        ]
        self.add_option(options)
        self.revise_infiles()
        self.set_options(self._sheet.options())

    def run(self):
        self.run_main()
        super(ReadsHumannWorkflow, self).run()

    def run_main(self):
        self.main = self.add_module("toolapps.humann_reads")
        if self.option("in_fastq"):
            qc = True
            options = {
                "in_fastq": self.option('in_fastq'),
                "qc": qc,
                "rm_host": self.option('rm_host'),
                "ref_database": self.option('ref_database'),
                "second_ref": self.option('second_ref'),
                "ref_undefined": self.option('ref_undefined'),
                "ref_undefined_name": self.option('ref_undefined_name'),
                "clean_fq": self.option('clean_fq'),
                "search_mode": self.option('search_mode'),
                "prescreen_threshold": self.option('prescreen_threshold'),
                "identity_threshold": self.option('identity_threshold'),
                "subject_coverage": self.option('subject_coverage'),
                "query_coverage": self.option('query_coverage'),
                "pathways": self.option('pathways'),
                "translated_alignment": self.option('aligin'),
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
        api_main.add_main_table("reads_humann2", main_id = self.option('main_id'))
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
            [".", "", "humann2结果输出目录",0],
        ])
        super(ReadsHumannWorkflow, self).end()