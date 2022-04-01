# -*- coding: utf-8 -*-
# __author__ = 'shaohua.yuan'
# last_modifiy = modified 2018.06.05

import os, re, subprocess, glob
from biocluster.agent import Agent
from biocluster.tool import Tool
import pandas as pd
import numpy as np
from biocluster.core.exceptions import OptionError
import traceback


class PathwayClusterAgent(Agent):
    """
    通路聚类分析
    """

    def __init__(self, parent):
        super(PathwayClusterAgent, self).__init__(parent)
        options = [
            {'name': 'exp', 'type': 'infile', 'format': 'metabolome.express,sequence.profile_table'},  # 表达矩阵文件
            {'name': 'sct', 'type': 'string', 'default': ''},  # 样本聚类算法，hierarchy、无
            {'name': 'scd', 'type': 'string', 'default': ''},  # 样本距离计算方式
            {'name': 'scm', 'type': 'string', 'default': ''},  # 样本聚类方式,"complete","average","single"
            {'name': 'mct', 'type': 'string', 'default': ''},  # 聚类算法，hierarchy、kmeans、无
            {'name': 'mcd', 'type': 'string', 'default': ''},  # 距离计算方式
            {'name': 'n_cluster', 'type': 'int', 'default': 0},  # 聚类数目
            {'name': 'mcm', 'type': 'string', 'default': ''},  # 聚类方式, hierarchy时使用，"complete","average","single"
        ]
        self.add_option(options)

    def check_options(self):
        """
        检查参数
        """
        if not self.option("exp").is_set:
            raise OptionError("请传入丰度矩阵！", code="34701301")
        if self.option("mct") == "kmeans":
            if not self.option("n_cluster") > 0:
                raise OptionError("使用kmeans时必须聚类数目必须大于1！", code="34701302")
        if self.option("sct") == "hierarchy":
            if not self.option("scm"):
                raise OptionError("样本使用hierarchy时必须传入层次聚类方式！", code="34701303")
        if self.option("mct") == "hierarchy":
            if not self.option("mcm"):
                raise OptionError("代谢物使用hierarchy时必须传入层次聚类方式！", code="34701304")
        '''
        if self.option("sct"):
            if not self.option("scd"):
                raise OptionError("必须输入样本距离计算方式！")
        if self.option("mct"):
            if not self.option("mcd"):
                raise OptionError("必须输入代谢物距离计算方式！")
        '''

    def set_resource(self):
        """
        所需资源
        """
        self._cpu = 3
        self._memory = '8G'

    def end(self):
        super(PathwayClusterAgent, self).end()


class PathwayClusterTool(Tool):
    """
    version 1.0
    """
    def __init__(self, config):
        super(PathwayClusterTool, self).__init__(config)
        software_dir = self.config.SOFTWARE_DIR
        self.python_path = 'program/Python/bin/python'
        self.script = self.config.PACKAGE_DIR + "/metabolome/scripts/corr_cluster.py"
        self.gcc = software_dir + '/gcc/5.1.0/bin'
        self.gcc_lib = software_dir + '/gcc/5.1.0/lib64'
        self.set_environ(PATH=self.gcc, LD_LIBRARY_PATH=self.gcc_lib)
        self.r_path = software_dir + "/program/R-3.3.1/bin:$PATH"
        self._r_home = software_dir + "/program/R-3.3.1/lib64/R/"
        self._LD_LIBRARY_PATH = software_dir + "/program/R-3.3.1/lib64/R/lib:$LD_LIBRARY_PATH"
        self.set_environ(PATH=self.r_path, R_HOME=self._r_home, LD_LIBRARY_PATH=self._LD_LIBRARY_PATH)

    def run(self):
        """
        运行
        """
        super(PathwayClusterTool, self).run()
        self.logger.info("开始运行命令！")
        self.run_cluster()
        self.set_output()

    def run_cluster(self):
        """
        """
        table = self.option("exp").prop["path"]
        cmd = '{} {} -exp {} -out {}'.format(self.python_path, self.script, table, self.work_dir)
        if self.option("sct"):
            cmd += " -sct " + self.option("sct")
        else:
            cmd += " --nsc"
        if self.option("mct"):
            cmd += " -gct " + self.option("mct")
        else:
            cmd += " --ngc"
        if self.option("scd"):
            cmd += " -scd " + self.option("scd")
        if self.option("mcd"):
            cmd += " -gcd " + self.option("mcd")
        cmd += " -n_clusters " + str(self.option("n_cluster"))
        if self.option("sct") == "hierarchy":
            cmd += " -scm " + self.option("scm")
        if self.option("mct") == "hierarchy":
            cmd += " -gcm " + self.option("mcm")
        command = self.add_command("cluster", cmd)
        command.run()
        self.wait()
        if command.return_code == 0:
            self.logger.info("cluster succeed")
        else:
            self.set_error("cluster failed!", code="34701301")
            raise Exception("cluster failed!")

    def set_output(self):
        self.logger.info("set output")
        if self.option("sct") == 'hierarchy':
            self.link_file("col.cluster_tree.txt", "metabset.cluster_tree.xls")
        if self.option("mct") == "hierarchy":  #kmean 时没树文件
            self.link_file("row.cluster_tree.txt", "pathway.cluster_tree.xls")
            self.link_file("row.cluster.txt", 'pathway.cluster.txt')
        else:
            self.link_file("row.kmeans_cluster.txt", 'pathway.cluster.txt')

        subclusters = glob.glob(self.work_dir + '/*subcluster*')
        for each in subclusters:
            each = each.split("/")[-1]
            newname = each.replace("row", "pathway")
            self.link_file(each, newname)

        self.end()

    def link_file(self, oldfile, newfile):
        oldfile = os.path.join(self.work_dir, oldfile)
        newfile = os.path.join(self.output_dir, newfile)
        if os.path.exists(newfile):
            os.remove(newfile)
        if os.path.exists(oldfile):
            os.link(oldfile, newfile)

