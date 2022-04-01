# !/mnt/ilustre/users/sanger-dev/app/program/Python/bin/python
# -*- coding: utf-8 -*-

from biocluster.agent import Agent
from biocluster.tool import Tool
from biocluster.core.exceptions import OptionError
from biocluster.config import Config
import os
import re
import shutil
import MySQLdb


class TableMergeAgent(Agent):
    """
    该tool用于将连接limus系统，下载客户信息表，并同步到monogo中，同时合并上机表
    version v1.0
    author: hongdongxuan
    last_modify: 20171013
    """
    def __init__(self, parent):
        super(TableMergekAgent, self).__init__(parent)
        options = [
            {"name": "table_up", "type": "infile", 'format': 'paternity_test.split_check'},  # 输入线上上机表
            {"name": "table_down", "type": "infile", "format": "paternity_test.split_check"},  # 输入线下上机表
            {"name": "batch_table", "type": 'infile', 'format': 'paternity_test.tab'}   # 输入批次表
        ]
        self.add_option(options)
        self.step.add_steps("TableMerge")
        self.on('start', self.stepstart)
        self.on('end', self.stepfinish)

    def stepstart(self):
        self.step.TableMerge.start()
        self.step.update()

    def stepfinish(self):
        self.step.TableMerge.finish()
        self.step.update()

    def check_options(self):
        """
        重写参数检测函数
        :return:
        """
        if not self.option("table_up").is_set:
            raise OptionError("必须输入线上上机表")
        if not self.option("table_down").is_set:
            raise OptionError("必须输入线下上机表")
        if not self.option('batch_table').is_set:
            raise OptionError("必须数据批次表！")
        return True

    def set_resource(self):
        """
        设置所需资源，需在之类中重写此方法 self._cpu ,self._memory
        :return:
        """
        self._cpu = 2
        self._memory = '10G'

    def end(self):
        super(TableMergeAgent, self).end()


class TableMergeTool(Tool):
    """
    运行获取客户信息，然后导入到mongo表中
    """
    def __init__(self, config):
        super(TableMergeTool, self).__init__(config)
        self._version = '1.0.1'

    def table_merge(self):
        """
        用于运行table合并以及样本名检查
        :return:
        """
        pass

    def pt_sample_name_check(self):
        """
        检查样本是否重命名以及检查样本名是否在批次表中有，以及样本名是否能够与客户信息表对应上
        :return:
        """
        self.api.paternity_test_v2.import_batch_table(self.option('batch_table').prop['path'])


    def run(self):
        super(TableMergeTool, self).run()
        self.pt_sample_name_check()
        self.table_merge()
        self.end()


class TableMerge(object):
    """
    用于合并线上与线下上机表
    """
    library = {"NIPT文库构建实验流程": "WS", "WQ gDNA建库实验流程_多重": "WQ", "WQ cfDNA建库实验流程_杂捕": "WQ",
               "FA gDNA建库实验流程_杂捕": "FA", "FA ctDNA建库实验流程_杂捕": "FA", "FA gDNA建库实验流程": "FA",
               "FA cfDNA建库实验流程_杂捕": "FA", "YCZL 建库实验流程": "YCZL", "QA 建库实验流程": "QA", "CT 建库实验流程": "CT",
               "HLY 建库实验流程": "HLY", "RXA 建库实验流程": "RXA", "snp 建库实验流程": "SNP", "XA RNA建库实验流程": "XA",
               "XA DNA建库实验流程": "XA", "甲基化 建库实验流程": "METH", "small RNA建库实验流程": "small"}

    analysistype = {"NIPT文库构建实验流程": "nipt", "WQ gDNA建库实验流程_多重": "dcpt", "WQ cfDNA建库实验流程_杂捕": "pt",
                    "FA gDNA建库实验流程_杂捕": "ctdna", "FA ctDNA建库实验流程_杂捕": "ctdna", "FA gDNA建库实验流程": "ctdna",
                    "FA cfDNA建库实验流程_杂捕": "ctdna", "YCZL 建库实验流程": "genetic_tumor", "QA 建库实验流程": "QA_str",
                    "CT 建库实验流程": "ct_str", "HLY 建库实验流程": "dc_str", "RXA 建库实验流程": "", "snp 建库实验流程": "",
                    "XA RNA建库实验流程": "blood_cancer", "XA DNA建库实验流程": "blood_cancer", "甲基化 建库实验流程": "",
                    "small RNA建库实验流程": ""}

    sampletype = {"全血": "QX", "血浆": "XJ", "蜡块": "SL", "石蜡": "SL", "石蜡切片": "SL", "石蜡切片(白片)": "SL",
                  "胸腹水": "XS", "手术标本": "XZ", "穿刺标本": "XZ", "穿刺样本": "XZ", "组织标本": "XZ", "新鲜组织": "XZ",
                  "蜡卷": "SL", "精斑": "JB", "亲子父本全血": "QQ", "指甲": "ZJ"}

    def __init__(self, up_table, down_table=None):
        self.up_table = up_table
        self.down_table = down_table
