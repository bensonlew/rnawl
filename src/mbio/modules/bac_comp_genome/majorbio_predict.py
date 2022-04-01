#!/usr/bin/env python
# -*- coding: utf-8 -*-
# __author__ == 'qingchen.zhang'

import os
import re
from biocluster.core.exceptions import OptionError
from biocluster.module import Module
import gevent
from biocluster.config import Config
import json
from collections import defaultdict
import time
import datetime
import shutil


class MajorbioPredictModule(Module):
    """
    微生物比较基因组从数据库中进行基因预测模块
    """
    def __init__(self, work_id):
        super(MajorbioPredictModule, self).__init__(work_id)
        options = [
            {"name": "sample", "type": "string"},###样本名称及对应的genome{"GCF_000217215.1": "NC_005823.1", "GCF_000217215.1": "NC_005867.1"}
            {"name": "genome_type", "type": "string"},###完成图基因组的类型，是chromosome还是plasmid，扫描图可以不同此字段；形式{"NC_005823.1": 'chromosome', "NC_005867.1": 'plasmid'}
            {"name": "assembly_type", "type": "string", "default": "draft"},###组装序列是完成图还是扫描图
        ]
        self.add_option(options)
        self.tools = []
        self.merge_tools = []

    def check_options(self):
        """
        参数检查
        :return:
        """
        self.logger.info("正在进行参数检查")
        if not self.option('assembly_type'):
            raise OptionError("必须提供样品类型")
        if not self.option('sample'):
            raise OptionError("必须提供样品名称")
        if not self.option('genome_type'):
            raise OptionError("必须提供基因组类型")

    def get_params(self):
        """
        获取参数
        :return:
        """
        sample = json.loads(self.option('sample'))
        self.genome_type = json.loads(self.option('genome_type'))
        self.genome_list = []
        sample_name = []
        for key in sample.keys():
            genome_id = sample[key]
            genome_list = genome_id.split(',')
            self.genome_list = sorted(genome_list)
            sample_name.append(key)
        name = ''.join(set(sample_name))
        return name

    def run_majorbio(self):
        """
        并行运行基因组预测
        :return:
        """
        self.sample = self.get_params()
        for genome in self.genome_list:
            self.majorbio = self.add_tool('bac_comp_genome.mongo_gene_predict')
            opts = {
                'sample': self.sample,
                'genome': genome,
                'genome_type': self.option('assembly_type'),
            }
            self.majorbio.set_options(opts)
            self.tools.append(self.majorbio)

    def merge_file(self):
        """
        将完成图的质粒和染色体所有序列合并起来
        :return:
        """
        if os.path.exists(self.work_dir + '/predict_dir'):
            shutil.rmtree(self.work_dir + '/predict_dir')
        os.mkdir(self.work_dir + '/predict_dir')
        for tool in self.tools:
            for file in os.listdir(tool.output_dir):
                old_path = os.path.join(tool.output_dir, file)
                new_path = os.path.join(self.work_dir + '/predict_dir', file)
                if os.path.exists(new_path):
                    os.remove(new_path)
                os.link(old_path, new_path)
        merge_dir = self.work_dir + '/predict_dir'
        merge_tool = self.add_tool('bac_comp_genome.cat_table')
        gff_list = ['_cds.gff', '_rrna.gff','_trna.gff',]
        seq_list = ['_cds.faa', '_cds.fnn', '_16s.fnn', '_rrna.fnn', '_trna.fnn']
        for file in gff_list:
            self.logger.info(file)
            opts = {
                'merge_dir': merge_dir,
                'prefix': self.sample,
                'header': 'True',
                'table_name': file
                }
            merge_tool.set_options(opts)
            self.merge_tools.append(merge_tool)
        for fnn in seq_list:
            self.logger.info(file)
            opts = {
                'merge_dir': merge_dir,
                'prefix': self.sample,
                'header': 'False',
                'table_name': fnn
                }
            merge_tool.set_options(opts)
            self.merge_tools.append(merge_tool)

    def predict_stat(self):
        """
        对所有的gff文件进行个数统计
        :return:
        """
        if len(self.genome_type.keys()) > 1:
            merge_dir = self.get_all_file()
        else:
            merge_dir = self.majorbio.output_dir



    def get_all_file(self):
        """
        获得前面所有merge后的结果文件
        :return:
        """
        merge_dir = self.work_dir + '/merge_dir'
        if os.path.exists(merge_dir):
            shutil.rmtree(merge_dir)
        os.mkdir(merge_dir)
        for tool in self.merge_tools:
            for file in os.listdir(tool.output_dir):
                old_path = os.path.join(tool.output_dir, file)
                new_file = os.path.join(merge_dir, file)
                if os.path.exists(new_file):
                    os.remove(new_file)
                os.link(old_path, new_file)
        return merge_dir

    def set_output(self):
        """
        设置结果文件目录
        :return:
        """

    def run(self):
        """
        运行
        :return:
        """
        super(MajorbioPredictModule, self).run()
        self.run_majorbio()
        if len(self.tools) == 1:
            self.on_rely(self.tools, self.predict_stat)
        elif len(self.tools) > 1:
            self.on_rely(self.tools, self.merge_file)
            if len(self.merge_tools) >= 1:
                self.on_rely(self.merge_tools, self.predict_stat)
            else:
                self.set_error('未能成功运行merge，请检查')
        for tool in self.tools:
            tool.run()
            gevent.sleep(0)

    def end(self):
        """
        结束
        :return:
        """
        super(MajorbioPredictModule, self).end()

