#!/usr/bin/env python
# -*- coding: utf-8 -*-
import os,re
import shutil
import glob
from biocluster.core.exceptions import OptionError
from biocluster.module import Module
from mbio.packages.bac_comp_genome.common_function import link_dir,link_file


class IslandModule(Module):
    """
    单个基因组基因组岛的预测
    last_modify: 2019.10.08
    """
    def __init__(self, work_id):
        super(IslandModule, self).__init__(work_id)
        options = [
            {"name": "fa_dir", "type": "infile", "format": "sequence.fasta_dir"},  #
            {"name": "gbk_dir", "type": "infile", "format": "gene_structure.gbk_dir"},  #
            {"name": "sample_name", "type": "string"},
            {"name": "gff", "type": "infile", "format": "gene_structure.gff3"},
            {"name": "genome", "type": "infile", "format": "sequence.fasta"},
        ]
        self.add_option(options)
        self.diomb = self.add_tool('bacgenome.island_dimob')
        self.islander = self.add_tool('bacgenome.island_islander')
        self.island = self.add_tool('bac_comp_genome.island_stat')
        self.list = [self.diomb, self.islander]

    def check_options(self):
        """
        检查参数
        """
        if not self.option('gbk_dir').is_set:
            raise OptionError("请设置基因组基因gbk文件夹不存在！")
        if not self.option('fa_dir').is_set:
            raise OptionError("请设置基因组学列文件夹！")

    def run_island_dimob(self):
        """
        island_dimob运行
        :return:
        """
        self.logger.info("正在island_dimob开始")
        self.diomb.set_options({
            'gbk_dir': self.option('gbk_dir')
        })
        self.diomb.run()

    def run_island_islander(self):
        """
        island_islander运行
        :return:
        """
        self.logger.info("正在island_islander开始")
        self.islander.set_options({
            'fa_dir': self.option("fa_dir")
        })
        self.islander.run()

    def run_island_stat(self):
        """
        对合并的island进行统计
        :return:
        """
        self.logger.info("正在对island数据统计开始")
        self.island.set_options({
            'diomb': self.diomb.option('out'),
            'islander': self.islander.option('out'),
            'gff': self.option('gff'),
            'sample_name': self.option("sample_name"),
            "genome": self.option("genome")
            })
        self.island.on("end", self.set_output)
        self.island.run()

    def run(self):
        super(IslandModule, self).run()
        self.run_island_dimob()
        self.run_island_islander()
        self.on_rely(self.list, self.run_island_stat)

    def set_output(self):
        self.logger.info("设置结果目录")
        if len(os.listdir(self.island.output_dir)) >=1:
            link_dir(self.island.output_dir, self.output_dir)
        self.end()

    def end(self):
        super(IslandModule, self).end()