# -*- coding: utf-8 -*-
# __author__ = 'shaohua.yuan'

from biocluster.agent import Agent
from biocluster.tool import Tool
import numpy as np
from biocluster.core.exceptions import OptionError
import os, re, glob
import subprocess
from mbio.packages.metabolome.common import Relation
import pandas as pd


class CorrTreeAgent(Agent):
    """
    计算相关系数和及相关性聚类树的工具
    version 1.0
    author: qindanhua
    last_modify: 2016.07.11
    """

    def __init__(self, parent):
        super(CorrTreeAgent, self).__init__(parent)
        options = [
            {'name': 'exp', 'type': 'infile', 'format': 'sequence.profile_table'},  # 表达矩阵文件
            {'name': 'sct', 'type': 'string', 'default': 'hierarchy'},  # 聚类算法，hierarchy、kmeans
            {'name': 'scd', 'type': 'string', 'default': 'euclidean'},  # 距离计算方式
            {'name': 'corr_method', 'type': 'string', 'default': 'pearson'},  # 相关性计算方法,"spearman", "pearson", "kendall"
            {'name': 'n_cluster', 'type': 'int', 'default': 0},  # 聚类数目，kmeans时使用
            {'name': 'scm', 'type': 'string', 'default': 'average'},  # 聚类方式, hierarchy时使用，"complete","average","single",""
            #{'name': 'metab_trans', 'type': 'infile', "format": "sequence.profile_table"},
            # 结果文件转化id使用，第一列为old_name,第二列为最后结果要使用的id, 可选参数
            {'name': 'file_tran', 'type': 'bool', 'default': True},
            {'name': 'pvalue', "type": "float", "default": 0.05}
        ]
        self.add_option(options)

    def check_options(self):
        """
        检查参数
        """
        if not self.option("exp").is_set:
            raise OptionError("请传入丰度矩阵！", code="32707001")
        if self.option("sct") == "kmeans":
            if not self.option("n_cluster") > 0:
                raise OptionError("使用kmeans时必须聚类数目必须大于1！", code="32707002")
        if self.option("sct") == "hierarchy":
            if not self.option("scm") or not self.option("scd"):
                raise OptionError("使用hierarchy时必须传入层次聚类方式和距离计算方式！", code="32707003")
        elif self.option("sct") == "kmeans":
            if not self.option("scd") or not self.option("n_cluster"):
                raise OptionError("使用kmeans时必须传入子聚类数目和距离计算方式！", code="32707004")

    def set_resource(self):
        """
        所需资源
        """
        self._cpu = 3
        self._memory = '5G'

    def end(self):
        super(CorrTreeAgent, self).end()


class CorrTreeTool(Tool):
    """
    version 1.0
    """

    def __init__(self, config):
        super(CorrTreeTool, self).__init__(config)
        software_dir = self.config.SOFTWARE_DIR
        self.python_path = 'miniconda2/bin/python'
        self.script = self.config.PACKAGE_DIR + "/meta/scripts/corr_cluster.py"
        self.gcc = software_dir + '/gcc/5.1.0/bin'
        self.gcc_lib = software_dir + '/gcc/5.1.0/lib64'
        self.set_environ(PATH=self.gcc, LD_LIBRARY_PATH=self.gcc_lib)
        self.r_path = software_dir + "/program/R-3.3.1/bin:$PATH"
        self._r_home = software_dir + "/program/R-3.3.1/lib64/R/"
        self._LD_LIBRARY_PATH = software_dir + "/program/R-3.3.1/lib64/R/lib:$LD_LIBRARY_PATH"
        self.set_environ(PATH=self.r_path, R_HOME=self._r_home, LD_LIBRARY_PATH=self._LD_LIBRARY_PATH)

    def correlation(self):
        table = self.option("exp").prop["path"]
        cmd = '{} {} -exp {} -out {} --corr --ngc '.format(self.python_path, self.script, self.table, self.work_dir)
        cmd += " -corr_method " + self.option("corr_method")

        if self.option("file_tran"):
            cmd += " --T "
        if self.option("sct"):
            cmd += " -sct " + self.option("sct")
        else:
            cmd += " --nsc "
        if self.option("sct") == "kmeans":
            cmd += " -n_cluster " + str(self.option("n_cluster"))
            cmd += " -scd " + self.option("scd")
        if self.option("sct") == "hierarchy":
            cmd += " -scm " + self.option("scm")
            cmd += " -scd " + self.option("scd")

        #cmd += ' -pvalue %s' % self.option('pvalue')
        command = self.add_command("correlation", cmd)
        command.run()
        self.wait()
        if command.return_code == 0:
            self.logger.info("correlation运行完成")
        else:
            self.set_error("correlation运行出错!", code="32707001")

    def set_output(self):
        self.logger.info("set output")

        if self.option("sct"):
            old_tree = os.path.join(self.work_dir, "corr.cluster_tree.txt")
            if os.path.exists(old_tree):
                self.link_file("corr.cluster_tree.txt", "corr.cluster_tree.xls")
        self.link_file("corr.xls", "corr.xls")
        self.link_file("pvalue.xls", "pvalue.xls")


    def link_file(self, oldfile, newfile):
        oldfile = os.path.join(self.work_dir, oldfile)
        newfile = os.path.join(self.output_dir, newfile)
        if os.path.exists(newfile):
            os.remove(newfile)
        os.link(oldfile, newfile)

    def check_otu_file(self):
        """
        对输入的表达矩阵文件进行检查，因为如果某个物种或者功能在各个样本的值相同的话，在用相关性系数方法计算时会报错， 那么这种情况就过滤掉，避免后面在计算时报错
        add by qingchen.zhang @20200116
        """
        table = self.option("exp").prop["path"]
        self.table = os.path.join(self.work_dir, "otu_table.xls")
        with open(table, 'r') as f, open(self.table, 'w') as w:
            lines = f.readlines()
            head = lines[0]
            w.write(head)
            if len(lines) <=2:
                self.set_error("The number of rows is less than 2, unable to calculate the correlation! ")
            for line in lines[1:]:
                line = line.strip().split("\t")
                sp_name = line[0]
                sp_abundance_list = line[1:]
                if len(set(sp_abundance_list)) == 1:##判断所有样本的丰度是否是相同的，如果相同就过滤掉
                    self.logger.info("{}: 物种在不同的样本中不能相同，需要过滤掉:{}".format(sp_name, sp_abundance_list))
                else:
                    w.write("\t".join(line) + "\n")

    def run(self):
        """
        运行
        """
        super(CorrTreeTool, self).run()
        self.check_otu_file()
        self.correlation()
        # self.correlation()
        # self.plot_hcluster()
        self.set_output()
        self.end()
