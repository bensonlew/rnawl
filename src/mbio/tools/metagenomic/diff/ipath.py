# -*- coding: utf-8 -*-
# __author__ = 'guhaidong'

from biocluster.agent import Agent
from biocluster.tool import Tool
import os
import re
from collections import defaultdict
import gridfs
from bson import ObjectId
import time
import json
import unittest
import pandas as pd
from mbio.packages.metagenomic.ipath import Ipath
from mbio.packages.metagenomic.common import styles
from biocluster.core.exceptions import OptionError
from biocluster.config import Config
import subprocess


class IpathAgent(Agent):
    """
    宏基因组IPATH 分析
    last_modify: 2018.10.18
    """

    def __init__(self, parent):
        super(IpathAgent, self).__init__(parent)
        options = [
            {"name": "ko_profile", "type": "infile", "format": "sequence.profile_table"},  # 基因丰度表
            # {"name": "anno_overview", "type": "infile", "format": "sequence.profile_table"},  # 基因集总览表
            {"name": "group_table", "type": "infile", "format": "meta.otu.group_table"},
            {"name": "color_style", "type": "string", "default": "default"}
        ]
        self.add_option(options)
        self.step.add_steps("ipath")
        self.on('start', self.stepstart)
        self.on('end', self.stepfinish)

    def stepstart(self):
        self.step.ipath.start()
        self.step.update()

    def stepfinish(self):
        self.step.ipath.finish()
        self.step.update()

    def check_options(self):
        """
        重写参数检测函数
        :return:
        """
        if not self.option("ko_profile").is_set:
            raise OptionError("请传入蛋白丰度表", code="32800301")
        if not self.option("group_table").is_set:
            raise OptionError("请传入分组表", code="32800302")
        try:
            styles[self.option('color_style')]
        except:
            raise OptionError("颜色类型设置不正确", code="32800303")

        return True

    def set_resource(self):
        """
        设置所需资源，需在之类中重写此方法 self._cpu ,self._memory
        :return:
        """
        self._cpu = 1
        self._memory = '5G'

    def end(self):
        super(IpathAgent, self).end()


class IpathTool(Tool):
    def __init__(self, config):
        super(IpathTool, self).__init__(config)
        # self.python = '/program/Python/bin/'
        # self.map_path = self.config.SOFTWARE_DIR + "/bioinfo/annotation/scripts/map5.r"
        self.ipath_db = self.config.SOFTWARE_DIR + "/database/IPATH"
        # self.image_magick = self.config.SOFTWARE_DIR + "/program/ImageMagick/bin/convert"
        self.ipath_input = self.work_dir + '/ipath_input.xls'
        # self.gene_ipath_input = self.work_dir + '/gene_ipath_input.xls'
        self.sets = []  # 第一个sets是Share
        self.group_color = {}

    def run(self):
        """
        运行
        :return:
        """
        super(IpathTool, self).run()
        # self.get_kegg_pics()
        self.generate_ipath_file()
        self.logger.info("ipath_db: %s" % self.ipath_db)
        self.logger.info("sets: %s" % self.sets)
        self.logger.info("group_color: %s" % self.group_color)
        self.logger.info("ipath_input: %s" % self.ipath_input)
        try:
            Ipath1 = Ipath()
            Ipath1.set_db(self.ipath_db)
            Ipath1.set_legend(self.sets)
            Ipath1.set_color_map(self.group_color)
            Ipath1.get_K_color_width(self.ipath_input)
            Ipath1.map_file()
        except Exception,e:
            self.logger.info(e)
            self.set_error("ipath运行错误: %s", variables=(e), code="32800301")
        self.logger.info("ipath运行完毕")
        self.set_output()
        self.end()

    def generate_ipath_file(self):
        data = pd.read_table(self.option("ko_profile").path, index_col=0)
        group = pd.read_table(self.option("group_table").path, index_col=0)
        group_header = group.columns[0]
        mapping = group.to_dict()[group_header]
        result = data.groupby(mapping, axis=1).sum()
        self.sets = ["Share"] + result.columns.tolist()
        if len(self.sets) == 2:
            self.logger.error("分组不能少于2个")
        colors = styles[self.option("color_style")]().get_colors(len(self.sets))
        group_color = dict(zip(self.sets, colors))
        table_data = result.apply(self.generate_ipath_data, axis=1, args=(group_color,))
        table_data.dropna(axis=0, how="any", inplace=True)
        self.logger.info(table_data.head())
        self.logger.info(colors[0])
        self.logger.info(colors[1:])
        share_data = table_data[table_data[0] == colors[0]]
        group_data = table_data[table_data[0].isin(colors[1:])]
        table_data = pd.concat([group_data, share_data])
        table_data.to_csv(self.ipath_input, sep="\t", header=False, index=True)
        self.group_color = group_color

    def generate_ipath_data(self, df, group_color):
        if df.prod() > 0:
            return pd.Series((group_color[self.sets[0]],"W15"), name=df.name)
        else:
            for group in self.sets[1:]:
                if df[group] >0 and df.drop([group]).prod() == 0:
                    return pd.Series((group_color[group],"W15"), name=df.name)
            return pd.Series((), name=df.name)  # 防止apply的第一条记录返回None,导致pd.apply自动判断返回的数据格式错误 @ ghd 20190307

    def set_output(self):
        all_files = ['ipath_input.xls', 'Metabolic_pathways.svg', 'Regulatory_pathways.svg',
        # all_files = ['Metabolic_pathways.svg', 'Regulatory_pathways.svg',
                     'Biosynthesis_of_secondary_metabolities.svg']
        for each in all_files:
            fname = os.path.basename(each)
            link = os.path.join(self.output_dir, fname)
            if os.path.exists(link):
                os.remove(link)
            os.link(each, link)

