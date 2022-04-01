# -*- coding: utf-8 -*-
# __author__ = 'zhaozhigang'

import os
from biocluster.workflow import Workflow
import datetime
import unittest
import types
import shutil
from bson.objectid import ObjectId


class ReadsBlastWorkflow(Workflow):
    """
    基于reads的humann2分析

    """

    def __init__(self, wsheet_object):
        self._sheet = wsheet_object
        super(ReadsBlastWorkflow, self).__init__(wsheet_object)
        options = [
            {'name': 'in_fastq', 'type': 'infile', 'format': 'sequence.fastq_dir'},  # 输入的fq文件夹
            {'name': 'clean_fq', 'type': 'infile', 'format': 'sequence.fastq_dir'},  # 输入的fq文件夹
            {'name': 'qc', 'type': 'string', 'default': "true"},  # true,false,fasta
            {'name': 'rm_host', 'type': 'bool', 'default': False},  # 是否需要去除宿主
            {'name': 'ref_database', 'type': 'string', 'default': ''},  # 宿主参考序列库中对应的物种名，eg：E.coli ,B.taurus
            {'name': 'second_ref', 'type': 'string', 'default': ''},  # 参考数据库二级下拉框数据
            {'name': 'ref_undefined', "type": 'infile', 'format': 'sequence.fasta_dir'},
            {'name': 'ref_undefined_name', 'type': 'string', 'default': 'undefined'},  # 自定义参考宿主名称，适应页面参数
            {'name': 'fasta_type', 'type': 'string', 'default': 'nucleotide'},  # nucleotide,protein
            {'name': 'fasta_dir', 'type': 'infile', 'format': 'toolapps.fasta_dir'},
            {'name': 'blast', 'type': 'string', 'default': 'blastn'},  # blastn,blastp,blastx
            {'name': 'database', 'type': 'string', 'default': 'NT'},  # NT,NR,Swiss-Prot,Custom
            {'name': 'ref_dir', 'type': 'infile', 'format': 'sequence.fasta_dir'},  # Custom时输入文件
            {'name': 'top_num', 'type': 'string', 'default': '1'},
            {'name': 'align_len', 'type': 'string', 'default': '50'},
            {'name': 'identity', 'type': 'string', 'default': '30'},
            {'name': 'evalue', 'type': 'string', 'default': '1e-5'},
            {"name": "main_id", "type": "string"},
            {"name": "update_info", "type": "string"},
        ]
        self.add_option(options)
        self.revise_infiles()
        self.set_options(self._sheet.options())
        print self.option('fasta_dir')

    def run(self):
        self.run_main()
        super(ReadsBlastWorkflow, self).run()

    def run_main(self):
        self.main = self.add_module("toolapps.meta_blast")
        if self.option('qc') in ['true', 'True']:
            options = {
                "in_fastq": self.option('in_fastq'),
                "qc": self.option('qc'),
                "rm_host": self.option('rm_host'),
                "ref_database": self.option('ref_database'),
                "second_ref": self.option('second_ref'),
                "ref_undefined": self.option('ref_undefined'),
                "ref_undefined_name": self.option('ref_undefined_name'),
                "blast": self.option('blast'),
                "database": self.option('database'),
                "ref_dir": self.option('ref_dir'),
                "top_num": self.option('top_num'),
                "align_len": self.option('align_len'),
                "identity": self.option('identity'),
                "evalue": self.option('evalue'),
            }
        elif self.option('qc') in ['false', 'False']:
            options = {
                "clean_fq": self.option('clean_fq'),
                "qc": self.option('qc'),
                "rm_host": self.option('rm_host'),
                "ref_database": self.option('ref_database'),
                "second_ref": self.option('second_ref'),
                "ref_undefined": self.option('ref_undefined'),
                "ref_undefined_name": self.option('ref_undefined_name'),
                "blast": self.option('blast'),
                "database": self.option('database'),
                "ref_dir": self.option('ref_dir'),
                "top_num": self.option('top_num'),
                "align_len": self.option('align_len'),
                "identity": self.option('identity'),
                "evalue": self.option('evalue'),
            }
            print options
        elif self.option('qc') in ['fasta', 'Fasta']:
            options = {
                "fasta_dir": self.option('fasta_dir'),
                "qc": self.option('qc'),
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
        api_main.add_main_table("reads_blast", main_id = self.option('main_id'))
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
            [".", "", "blast分析结果文件目录",0],
        ])
        super(ReadsBlastWorkflow, self).end()