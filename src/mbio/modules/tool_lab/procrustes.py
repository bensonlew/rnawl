# -*- coding: utf-8 -*-
# __author__ = 'zhangyitong'

import os
import shutil
import unittest

import pandas as pd
import glob
from biocluster.core.exceptions import OptionError
from biocluster.module import Module
from mbio.files.sequence.file_sample import FileSampleFile


class ProcrustesModule(Module):
    """
    """

    def __init__(self, work_id):
        super(ProcrustesModule, self).__init__(work_id)
        options = [
            {'name': 'ref_table', 'type': 'infile', 'format': 'sequence.profile_table', 'required': True},
            {'name': 'query_table', 'type': 'infile', 'format': 'sequence.profile_table', 'required': True},
            {'name': 'group_file', 'type': 'infile', 'format': 'meta.otu.group_table'},
            {'name': 'dist', 'type': 'string', 'default': 'euclidean'},
            {'name': 'method', 'type': 'string'},
        ]
        self.add_option(options)
        # if self.option('method') == 'pca':
        #     self.pca_ref_tool = self.add_tool("tool_lab.procrustes.pca")
        #     self.pca_query_tool = self.add_tool("tool_lab.procrustes.pca")
        # elif self.option('method') == 'pcoa':
        #     self.discalc_ref_tool = self.add_tool("tool_lab.procrustes.distance_calc")
        #     self.discalc_query_tool = self.add_tool("tool_lab.procrustes.distance_calc")
        #     self.pcoa_ref_tool = self.add_tool("tool_lab.procrustes.pcoa")
        #     self.pcoa_query_tool = self.add_tool("tool_lab.procrustes.pcoa")
        self.procrustes_tool = self.add_tool("tool_lab.procrustes.procrustes")
        self.run_tools = []  # 并行运行的tool

    def check_options(self):
        """
        检查参数
        """
        if self.option('method') not in ["pca", "pcoa"]:
            raise OptionError("PARAMETERS ERROR:wrong value of method (%s), pca or pcoa expected!" % self.option('method'))
        return True

    def discalc_run_ref(self):
        """
        运行计算reference距离矩阵，代谢物丰度文件
        :return:
        """
        options = {
            "otutable": self.option('ref_table').prop['path'],
            "method": self.option("dist")
        }
        self.discalc_ref_tool.set_options(options)
        self.discalc_ref_tool.run()

    def discalc_run_query(self):
        """
        运行计算query距离矩阵，比如自定义上传的文件
        :return:
        """
        options = {
            "otutable": self.option('query_table').prop['path'],
            "method": self.option("dist")
        }
        self.discalc_query_tool.set_options(options)
        self.discalc_query_tool.run()

    def pcoa_run_ref(self):
        options = {
            'dis_matrix': self.discalc_ref_tool.option('dis_matrix'),
        }
        self.pcoa_ref_tool.set_options(options)
        self.pcoa_ref_tool.run()

    def pcoa_run_query(self):
        options = {
            'dis_matrix': self.discalc_query_tool.option('dis_matrix'),
        }
        self.pcoa_query_tool.set_options(options)
        self.pcoa_query_tool.run()

    def pca_run_ref(self):
        options = {
            'otutable': self.option('ref_table').prop['path'],
            'group_table': self.option('group_file').prop['path'],
        }
        self.pca_ref_tool.set_options(options)
        self.pca_ref_tool.run()

    def pca_run_query(self):
        options = {
            'otutable': self.option('query_table').prop['path'],
            'group_table': self.option('group_file').prop['path'],
        }
        self.pca_query_tool.set_options(options)
        self.pca_query_tool.run()

    def procrustes_run(self):
        coord_ref = ''
        coord_query = ''
        if self.option('method') == 'pca':
            coord_ref = os.path.join(self.pca_ref_tool.work_dir, "pca/pca.txt")
            coord_query = os.path.join(self.pca_query_tool.work_dir, "pca/pca.txt")
        elif self.option('method') == 'pcoa':
            coord_ref = os.path.join(self.pcoa_ref_tool.work_dir, "pcoa/pcoa.txt")
            coord_query = os.path.join(self.pcoa_query_tool.work_dir, "pcoa/pcoa.txt")

        options = {
            'coord_ref': coord_ref,
            'coord_query': coord_query
        }
        self.procrustes_tool.set_options(options)
        self.procrustes_tool.run()

    def set_output(self):
        target_files = glob.glob(self.procrustes_tool.output_dir + '/*transformed_reference.txt')
        target_files += glob.glob(self.procrustes_tool.output_dir + '/*transformed_q*.txt')
        target_files += glob.glob(self.procrustes_tool.output_dir + '/procrustes_results.txt')
        for each in target_files:
            name = os.path.basename(each)
            link = os.path.join(self.output_dir, name)
            if os.path.exists(link):
                os.remove(link)
            os.link(each, link)
        self.end()

    def run(self):
        super(ProcrustesModule, self).run()

        '''
        定义依赖关系和执行顺序
        '''
        if self.option('method') == 'pca':
            self.pca_ref_tool = self.add_tool("tool_lab.procrustes.pca")
            self.pca_query_tool = self.add_tool("tool_lab.procrustes.pca")
        elif self.option('method') == 'pcoa':
            self.discalc_ref_tool = self.add_tool("tool_lab.procrustes.distance_calc")
            self.discalc_query_tool = self.add_tool("tool_lab.procrustes.distance_calc")
            self.pcoa_ref_tool = self.add_tool("tool_lab.procrustes.pcoa")
            self.pcoa_query_tool = self.add_tool("tool_lab.procrustes.pcoa")
        if self.option('method') == 'pca':
            self.run_tools.append(self.pca_ref_tool)
            self.run_tools.append(self.pca_query_tool)
        elif self.option('method') == 'pcoa':
            self.run_tools.append(self.pcoa_ref_tool)
            self.run_tools.append(self.pcoa_query_tool)
            self.discalc_ref_tool.on('end', self.pcoa_run_ref)
            self.discalc_query_tool.on('end', self.pcoa_run_query)

        if self.option('method') == 'pca':
            self.pca_run_ref()
            self.pca_run_query()
        elif self.option('method') == 'pcoa':
            self.discalc_run_ref()
            self.discalc_run_query()
        self.on_rely(self.run_tools, self.procrustes_run)
        self.procrustes_tool.on("end", self.set_output)

    def end(self):
        super(ProcrustesModule, self).end()
