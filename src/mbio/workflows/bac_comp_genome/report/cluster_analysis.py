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



class ClusterAnalysisWorkflow(Workflow):
    def __init__(self, wsheet_object):
        """
        细菌比较基因组泛基因组分析的第一步：
        聚类、泛基因组大小计算和变异分析；
        第二步 对聚类结果进行分组方案的合并
        第三步 对聚类的结果进行venn图的分析
        :return:
        """
        self._sheet = wsheet_object
        super(ClusterAnalysisWorkflow, self).__init__(wsheet_object)
        options = [
            {"name": "analysis", "type": "string",
             "default": "pca,pcoa,nmds,hcluster,plsda"},
            {"name": "dis_method", "type": "string", "default": "bray_curtis"},
            # 当设定此值时，dbrda的计算方式将会改变，使用R中自带的距离算法，而不是先计算好距离矩阵，此处的计算方式与一般的距离计算的的值不一致
            {"name": "otutable", "type": "infile", "format": "meta.otu.otu_table, meta.otu.tax_summary_dir"},
            {"name": "level", "type": "string", "default": "otu"},
            {"name": "phy_newick", "type": "infile", "format": "meta.beta_diversity.newick_tree"},
            {"name": "permutations", "type": "int", "default": 999},
            {"name": "linkage", "type": "string", "default": "average"},
            {"name": "group", "type": "infile", "format": "meta.otu.group_table"},
            {"name": "grouplab", "type": "string", "default": ""},
            {"name": "dis_matrix", "type": "outfile", "format": "meta.beta_diversity.distance_matrix"},
        ]

        self.pca = self.add_module("bac_comp_genome.beta_diversity")
        self.heatmap = self.add_module("bac_comp_genome.hcluster_heatmap")
        self.IMPORT_REPORT_DATA = True
        self.IMPORT_REPORT_AFTER_END = False
        self.all_tools = []
        self.add_option(options)
        self.set_options(self._sheet.options())
        self.sample = []

    def check_options(self):
        """
        进行参数二次检查
        :return:
        """
        self.logger.info("开始pan的参数检查")
        return True



    def set_output(self):
        """
        设置结果文件目录
        :return:
        """
        self.logger.info('设置结果目录')

        self.logger.info('设置结果目录成功')
        self.set_db()

    def get_sample_list(self):
        """
        获取样本的list，便于导表，交互分析不用，直接传过来，工作流要写
        :return:
        """
        samples_li = []
        sample_dir = self.option("fasta_dir").prop['path']
        sample_list = os.listdir(sample_dir)
        for file in sample_list:
            if file.endswith('faa'):
                sample_name = file.strip(".faa")
                samples_li.append(sample_name)
        return samples_li

    def set_db(self):
        """
        将数据导入mongo
        :return:
        """
        self.logger.info('正在写入mongo数据库')


        self.end()


    def run(self):
        """
        开始运行了
        :return:
        """
        self.logger.info("开始运行")
        self.down_file()
        self.cluster.on('end', self.run_all_tool)
        self.run_cluster()
        super(ClusterAnalysisWorkflow, self).run()

    def end(self):
        """
        结束了
        :return:
        """
        super(ClusterAnalysisWorkflow, self).end()
