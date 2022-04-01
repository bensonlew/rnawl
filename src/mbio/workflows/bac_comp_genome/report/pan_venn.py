# -*- coding: utf-8 -*-
#__author__ = 'qingchen.zhang'
import os,shutil
import re
import types
from biocluster.workflow import Workflow
from bson import ObjectId
from biocluster.core.exceptions import OptionError
import datetime
from biocluster.file import download, exists
from collections import defaultdict
import shutil
import gevent
from mbio.packages.bac_comp_genome.common_function import link_dir



class PanVennWorkflow(Workflow):
    def __init__(self, wsheet_object):
        """
        细菌比较基因组泛基因组分析
        第三步 对聚类的结果进行venn图的分析
        :return:
        """
        self._sheet = wsheet_object
        super(PanVennWorkflow, self).__init__(wsheet_object)
        options = [
            {"name": "group_table", "type": "infile", "format": "meta.otu.group_table"},#输入group表选择哪些样本进行导表和计算
            {"name": "cluster_file", "type": "infile", "format": "sequence.profile_table"},  # 物种注释总览表
            {"name": "main_id", "type": "string"},#输入pan_venn表的id
            {"name": "update_info", "type": "string"},
            {'name': 'group_detail', 'type': 'string'}, ##to_file专用
        ]
        self.venn = self.add_tool("bac_comp_genome.pan_venn")
        self.IMPORT_REPORT_DATA = True
        self.IMPORT_REPORT_AFTER_END = False
        self.add_option(options)
        self.set_options(self._sheet.options())

    def check_options(self):
        """
        进行参数二次检查
        :return:
        """
        self.logger.info("开始pan的参数检查")
        if not self.option("group_table").is_set:
            raise OptionError("必须设置输入的group_table")
        if not self.option("cluster_file").is_set:
            raise OptionError("请提供输入序列的文件夹cluster_file！")
        return True

    def run_venn(self):
        """
        Venn图的分析
        :return:
        """
        self.logger.info("开始进行venn图分析")
        self.venn.set_options({
            "cluster": self.option("cluster_file"),
            "group_table": self.option("group_table"),
            "version": "v1"
            })
        self.venn.run()

    def set_output(self):
        """
        设置结果文件目录
        :return:
        """
        self.logger.info('设置结果目录')
        result_path = os.path.join(self.output_dir)
        link_dir(self.venn.output_dir, result_path)
        self.logger.info('设置结果目录成功')
        self.set_db()

    def set_db(self):
        """
        将数据导入mongo
        :return:
        """
        self.logger.info('正在写入mongo数据库')
        self.remote_dir = self._sheet.output
        api_venn = self.api.api('bac_comp_genome.pan_venn')
        self.logger.info("正在进行venn图分析")
        venn_path = os.path.join(self.output_dir, "venn_table.xls")
        venn_id = self.option("main_id")
        api_venn.add_venn_detail(venn_id, venn_path)
        self.end()

    def run(self):
        """
        开始运行了
        :return:
        """
        self.logger.info("开始运行")

        self.venn.on('end', self.set_output)
        self.run_venn()
        super(PanVennWorkflow, self).run()

    def end(self):
        """
        结束了
        :return:
        """
        self.logger.info("开始结果文件上传")
        sdir = self.add_upload_dir(self.output_dir)
        repaths = [
            ["PanVenn", "", "PanVenn分析结果输出目录",0,""],
            ["venn_table.xls", "xls", "PanVenn分析结果表",0,""],
            ]
        sdir.add_relpath_rules(repaths)
        super(PanVennWorkflow, self).end()
