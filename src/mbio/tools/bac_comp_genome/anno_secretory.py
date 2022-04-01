# -*- coding: utf-8 -*-
# __author__ = 'gaohao'
# last_modify : 2019.09.17

from biocluster.agent import Agent
from biocluster.tool import Tool
from biocluster.core.exceptions import OptionError
from mbio.packages.align.blast.xml2table import xml2table
import pandas as pd
import subprocess, os,re


class AnnoSecretoryAgent(Agent):
    """
    调用kegg注释需要的脚本 kegg_xlm_mongo.py 获取注释信息以及对ko_name进行取名规则处理
    """

    def __init__(self, parent):
        super(AnnoSecretoryAgent, self).__init__(parent)
        options = [
            {"name": "kegg_anno", "type": "infile", "format": "sequence.profile_table"},  # 输入的比对结果xml文件
        ]
        self.add_option(options)

    def check_options(self):
        if not self.option("kegg_anno").is_set:
            raise OptionError("必须设置kegg注释输入文件")
        return True

    def set_resource(self):
        self._cpu = 2
        self._memory = '5G'

    def end(self):
        super(AnnoSecretoryAgent, self).end()


class AnnoSecretoryTool(Tool):
    def __init__(self, config):
        super(AnnoSecretoryTool, self).__init__(config)
        self._version = "1.0"
        self.python_path = self.config.SOFTWARE_DIR + "/program/Python/bin/python"
        self.python_script = self.config.PACKAGE_DIR + '/bac_comp_genome/secretory_anno.py'
        self.ref = self.config.SOFTWARE_DIR + "/database/Secretory/Secretion_System.xls"

    def run(self):
        """
        运行
        :return:
        """
        super(AnnoSecretoryTool, self).run()
        self.run_secretory_anno()
        self.set_output()
        self.end()

    def run_secretory_anno(self):
        table = self.option('kegg_anno').prop['path']
        cmd = '{} {} -s {} -i {} -o {}'.format(self.python_path, self.python_script, self.ref, table, self.output_dir+ "/secretory_system.xls")
        self.logger.info(cmd)
        try:
            subprocess.check_output(cmd, shell=True)
            self.logger.info('运行run_secretory_anno_anno完成')
        except subprocess.CalledProcessError:
            self.set_error('运行run_secretory_anno_anno出错')

    def set_output(self):
        pass