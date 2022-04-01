# -*- coding: utf-8 -*-
# __author__ = 'shaohua.yuan'
# last_modifiy = modified 2018.06.05

import os, re, subprocess, glob
from biocluster.agent import Agent
from biocluster.tool import Tool
import pandas as pd
import numpy as np
from biocluster.core.exceptions import OptionError
from mbio.packages.metabolome.scripts.mul_diff_stat import MulDiffStat
import traceback
from mbio.packages.metabolome.common import Relation

class ExpPcaMetaAgent(Agent):
    """
    两组差异检验，student T检验，welch T检验，wilcox秩和检验
    """

    def __init__(self, parent):
        super(ExpPcaMetaAgent, self).__init__(parent)
        options = [
            {'name': 'exp_file', 'type': 'infile', 'format': 'labelfree.express_matrix'},  # 表达矩阵文件
            {"name": "group_file", "type": "infile", "format": "labelfree.common"},  # 分组文件
            {'name': 'mul_type', 'type': 'string', 'default': 'pca'},  # 多元统计类型，pca，plsda, oplsda ，可以是“pca;plsda”形式
            {'name': 'confidence', 'type': 'string', 'default': '0.95'},  # 置信度，与mul_type对应的“0.95;0.95”形式
            {'name': 'perm', 'type': 'string', 'default': ''},  # 置换次数，与mul_type对应的“200;100”形式
            {'name': 'data_trans', 'type': 'string', 'default': 'UV'},  # 数据转化方法："UV","Ctr","Par"，"", 与mul_type对应个数
        ]
        self.add_option(options)

    def check_options(self):
        """
        检查参数
        """
        if not self.option("exp_file").is_set:
            raise OptionError("请传入丰度矩阵！", code="34700101")
        if not self.option("group_file").is_set:
            raise OptionError("请传入group文件！", code="34700102")
        return True

    def set_resource(self):
        """
        所需资源
        """
        self._cpu = 3
        self._memory = '5G'

    def end(self):
        super(ExpPcaMetaAgent, self).end()


class ExpPcaMetaTool(Tool):
    """
    version 1.0
    """

    def __init__(self, config):
        super(ExpPcaMetaTool, self).__init__(config)
        #self.r_path = self.config.SOFTWARE_DIR + "/program/R-3.3.3/bin/Rscript"
        self.set_environ(R_LIBS=self.config.SOFTWARE_DIR + '/program/R-3.3.1/lib64/R/library')
        self.r_path = "/program/R-3.3.1/bin/Rscript"

    def run(self):
        """
        运行
        """
        super(ExpPcaMetaTool, self).run()
        self.logger.info("开始运行命令！")
        self.logger.info(self.option("exp_file").prop["path"])
        self.logger.info(self.option("group_file").prop["path"])
        self.run_mul_stat()
        self.set_output()


    def run_mul_stat(self):
        """
        pca、plsda、oplsda分析
        """
        diffstat = MulDiffStat()
        diffstat.mul_pca(self.option("exp_file").prop["path"], self.output_dir, self.option("group_file").prop["path"], ci=self.option("confidence"), data_trans=self.option("data_trans"), mul_type=self.option("mul_type"), perm=self.option("mul_type"))
        rname = glob.glob(os.path.join(self.work_dir, 'run_*.r'))[0]
        cmd = self.r_path + ' ' + rname
        command = self.add_command("pca", cmd).run()
        self.wait(command)
        if command.return_code == 0:
            self.logger.info("run {} succeed".format(cmd))
        else:
            self.set_error("run %s failed", variables=(cmd), code="34700101")
            raise Exception("run {} failed".format(cmd))


    def set_output(self):
        """
        id 转化
        """
        raw_dir = os.getcwd()
        os.chdir(self.output_dir)
        os.remove('PCA.ellipse.xls')
        with open('PCA.model.xls') as mr, open('pca_importance.xls', 'w') as iw:
            _ = mr.readline()
            iw.write('Sample_ID\tProportion of Variance\n')
            for line in mr:
                if not line.strip():
                    continue
                l = line.strip().split('\t')
                iw.write('%s\t%s\n' % (l[0].replace('p', 'PC'), l[1]))
        os.remove('PCA.model.xls')
        with open('PCA.sites.xls') as sr, open('pca_sites.xls', 'w') as sw:
            header = sr.readline()
            sw.write('Sample_ID' + header.replace('p', 'PC'))
            for line in sr:
                sw.write(line)
        os.remove('PCA.sites.xls')
        os.chdir(raw_dir)
        self.end()
