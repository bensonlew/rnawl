# -*- coding: utf-8 -*-
# __author__ = 'gaohao'

import pandas as pd
import os
from biocluster.workflow import Workflow
import types
from bson.objectid import ObjectId
from mbio.packages.tool_lab.common import down_seq_files
import shutil
import json
import HTMLParser


class PmlstWorkflow(Workflow):
    """
    单个基因组的PMLST分析
    """
    def __init__(self, wsheet_object):
        self.samples = {}
        self._sheet = wsheet_object
        if self._sheet.option('project_data'):
            self.project_data = eval(HTMLParser.HTMLParser().unescape(self._sheet.option('project_data')))
            for i in self.project_data['specimens']:
                self.samples[i['id']] = i['name']
            assemble_dir = os.path.join(self._sheet.work_dir, "assemble_dir")
            if os.path.exists(assemble_dir):
                shutil.rmtree(assemble_dir)
            (self.assemble_dir, self.analysis_type) = down_seq_files(self.project_data['my_type'],
                                                                     self.project_data['db_version'], assemble_dir,
                                                                     self._sheet.option("project_task_id"),
                                                                     self.samples)
        super(PmlstWorkflow, self).__init__(wsheet_object)
        options = [
            {"name": "fasta", "type": "infile", "format": "sequence.fasta"},# 序列文件
            {"name": "species", "type": "string", },# 物种名称
            {'name': 'task_id', 'type': 'string'},
            {'name': 'project_data','type':'string'},
            {'name': 'project_task_id', 'type': 'string'},
            {'name': 'update_info', 'type': 'string'},
            {'name': 'main_id', 'type': 'string'}
        ]
        self.add_option(options)
        self.revise_infiles()
        self.set_options(self._sheet.options())

    def run(self):
        self.run_pmlst()
        super(PmlstWorkflow, self).run()

    def download_file(self):
        """
        download file from s3
        :return:
        """
        if self.analysis_type in ['complete']:
            assemble_dir = os.path.join(self.work_dir, "assemble_dir")
            for sample in self.samples.keys():
                sample_path = os.path.join(assemble_dir, self.samples[sample] + ".fna")
                assemble_dir3 = os.path.join(assemble_dir, self.samples[sample])
                dir_list = os.listdir(assemble_dir3)
                for file2 in dir_list:
                    os.system("cat {} >> {}".format(os.path.join(assemble_dir3, file2), sample_path))

    def run_pmlst(self):
        self.pmlst = self.add_tool("tool_lab.pmlst")
        if self.option('fasta').is_set:
            sample = ".".join(os.path.basename(self.option('fasta').prop['path']).split(".")[0:-1])
            self.pmlst.set_options({
                'fasta': self.option('fasta'),
                "sample_name": sample,
                'species': self.option('species')
            })
        if self.option('project_data'):
            self.download_file()
            sample = self.samples.values()[0]
            fasta = os.path.join(self.work_dir, "assemble_dir", sample + ".fna")
            self.pmlst.set_options({
                'fasta': fasta,
                "sample_name": sample,
                'species': self.option('species')
            })
        self.pmlst.on('end', self.set_db)
        self.pmlst.run()


    def set_db(self):
        """
        导表
        """
        sample = ''
        if self.option('fasta').is_set:
            sample = ".".join(os.path.basename(self.option('fasta').prop['path']).split(".")[0:-1])
        if self.option('project_data'):
            sample = self.samples.values()[0]
        pmlst = self.api.api('tool_lab.pmlst')
        path1 =self.pmlst.output_dir + "/" + sample + ".pmlst.ST.xls"
        path2 = self.pmlst.output_dir + "/" + sample + ".pmlst.detail.xls"
        pmlst.add_mlst_detail(ObjectId(self.option("main_id")), path1, path2, sample)
        self.end()

    def end(self):
        result_dir = self.add_upload_dir(self.pmlst.output_dir)
        result_dir.add_relpath_rules([
            [".", "", "Pmlst结果输出目录"],
            ["./*.pmlst.detail.xls", "xls", "PMLST分型比对详情表"],
            ["./*.HitMLST.fasta", "xls", "PMLST分型比对等位基因信息"],
            ["./*.pmlst.ST.xls", "xls", "PMLST分型统计表"],
        ])
        result_dir.add_regexp_rules([
            ["./*.pmlst.detail.xls", "xls", "PMLST分型比对详情表"],
            ["./*.HitMLST.fasta", "xls", "PMLST分型比对等位基因信息"],
            ["./*.pmlst.ST.xls", "xls", "PMLST分型统计表"],
        ])
        super(PmlstWorkflow, self).end()
