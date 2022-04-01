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


class SerotypefinderWorkflow(Workflow):
    """
    多个基因组的Serotypefinder分析
    """
    def __init__(self, wsheet_object):
        self._sheet = wsheet_object
        super(SerotypefinderWorkflow, self).__init__(wsheet_object)
        options = [
            {"name": "fasta_dir", "type": "infile", "format": "tool_lab.fasta_dir"},# 序列文件夹
            {"name": "fasta", "type": "infile", "format": "sequence.fasta"},  # 序列
            {'name': 'task_id', 'type': 'string'},
            {'name': 'update_info', 'type': 'string'},
            {'name': 'type', 'type': 'string'},  ##是页面单个基因组还是多个基因组
            {'name': 'main_id', 'type': 'string'}
        ]
        self.add_option(options)
        self.revise_infiles()
        self.set_options(self._sheet.options())

    def run(self):
        self.run_kmerfinder()
        super(SerotypefinderWorkflow, self).run()

    def run_kmerfinder(self):
        if self.option('fasta').is_set:
            fasta_dir = self.get_dir(self.option('fasta').prop['path'])
        elif self.option('fasta_dir').is_set:
            fasta_dir = self.option('fasta_dir').prop['path']
        self.serotypefinder = self.add_module("tool_lab.serotypefinder")
        self.serotypefinder.set_options({
            'fasta_dir': fasta_dir,
        })
        self.serotypefinder.on('end', self.set_db)
        self.serotypefinder.run()

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
        link_dir(self.serotypefinder.output_dir, self.output_dir)
        if self.option('fasta').is_set:
            list_file = self.work_dir+"/seq/list.txt"
        elif self.option('fasta_dir').is_set:
            list_file = self.option('fasta_dir').prop['path']+"/list.txt"
        samples =self.get_sample(list_file)
        serotypefinder = self.api.api('tool_lab.serotypefinder')
        if len(os.listdir(self.serotypefinder.output_dir)) >0:
            for sample in samples:
                if os.path.exists(self.output_dir+"/"+sample+".SerotypeDetail.xls"):
                    serotypefinder.add_serotypefinder_detail(ObjectId(self.option("main_id")), self.output_dir+"/"+sample+".SerotypeDetail.xls", sample)
                if os.path.exists(self.output_dir+"/"+sample+".SerotypeStat.xls"):
                    serotypefinder.add_serotypefinder_stat(ObjectId(self.option("main_id")), self.output_dir+"/"+sample+".SerotypeStat.xls")
        self.end()


    def end(self):
        result_dir = self.add_upload_dir(self.output_dir)
        result_dir.add_relpath_rules([
            [".", "", "Mlst结果输出目录"],
            ["./*.SerotypeDetail.xls", "xls", "Serotype对详情表"],
            ["./*.HitMLST.fasta", "xls", "Serotype对等位基因信息"],
            ["./*.SerotypeStat.xls", "xls", "Serotype统计表"],
        ])
        result_dir.add_regexp_rules([
            ["./*.SerotypeDetail.xls", "xls", "Serotype比对详情表"],
            ["./*.HitMLST.fasta", "xls", "Serotype比对等位基因信息"],
            ["./*.SerotypeStat.xls", "xls", "Serotype统计表"],
        ])
        super(SerotypefinderWorkflow, self).end()

    def get_sample(self,file):
        list = []
        with open (file, "r") as f:
            lines = f.readlines()
            for line in lines[1:]:
                lin = line.strip().split("\t")
                list.append(lin[0])
        return list
