# -*- coding: utf-8 -*-
# __author__ = 'hao.gao'

import pandas as pd
import os
from biocluster.workflow import Workflow
import datetime
import types
from bson.objectid import ObjectId
import re
import json
import shutil
from mbio.packages.bac_comp_genome.common_function import link_dir,link_file


class KmerfinderWorkflow(Workflow):
    """
    多个基因组的kmerfinder分析
    """
    def __init__(self, wsheet_object):
        self._sheet = wsheet_object
        super(KmerfinderWorkflow, self).__init__(wsheet_object)
        options = [
            {"name": "fasta_dir", "type": "infile", "format": "tool_lab.fasta_dir"},# 序列文件夹
            {"name": "fasta", "type": "infile", "format": "sequence.fasta"},  # 序列
            {"name": "species", "type": "string", },# 物种名称
            {'name': 'task_id', 'type': 'string'},
            {'name': 'type', 'type': 'string'},  ##是页面单个基因组还是多个基因组
            {'name': 'update_info', 'type': 'string'},
            {'name': 'main_id', 'type': 'string'}
        ]
        self.add_option(options)
        self.revise_infiles()
        self.set_options(self._sheet.options())

    def run(self):
        self.run_kmerfinder()
        super(KmerfinderWorkflow, self).run()

    def run_kmerfinder(self):
        if self.option('fasta').is_set:
            fasta_dir = self.get_dir(self.option('fasta').prop['path'])
        elif self.option('fasta_dir').is_set:
            fasta_dir = self.option('fasta_dir').prop['path']
        self.kmerfinder = self.add_module("tool_lab.kmerfinder")
        self.kmerfinder.set_options({
            'fasta_dir': fasta_dir,
            'species': self.option('species')
        })
        self.kmerfinder.on('end', self.set_db)
        self.kmerfinder.run()

    def get_dir(self, file):
        sample = ".".join(os.path.basename(file).split('.')[0:-1])
        if os.path.exists(self.work_dir+"/seq"):
            shutil.rmtree(self.work_dir+"/seq")
        os.mkdir(self.work_dir+"/seq")
        self.logger.info(file)
        os.link(file, self.work_dir+"/seq/"+sample+".fasta")
        with open(self.work_dir+"/seq/list.txt", "w") as g:
            g.write("Sample\tfile\n")
            g.write("{}\t{}\n".format(sample, sample+".fasta"))
        return self.work_dir+"/seq/"

    def set_db(self):

        """
         保存结果标准化数据到mongo数据库中
        """
        kmerfinder = self.api.api('tool_lab.kmerfinder')
        if len(os.listdir(self.kmerfinder.output_dir)) >0:
            link_dir(self.kmerfinder.output_dir, self.output_dir)
            path1 = self.output_dir
            kmerfinder.add_kmerfinder_detail(ObjectId(self.option("main_id")), path1)
        self.end()


    def end(self):
        result_dir = self.add_upload_dir(self.output_dir)
        result_dir.add_relpath_rules([
            [".", "", "结果输出目录"],
            ["./*.KmerFinderResult.xls", "xls", "KmerFinder详情表"],
        ])
        result_dir.add_regexp_rules([
            ["./*.KmerFinderResult.xls", "xls", "KmerFinder详情表"],
        ])
        super(KmerfinderWorkflow, self).end()
