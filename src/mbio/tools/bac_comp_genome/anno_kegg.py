# -*- coding: utf-8 -*-
# __author__ = 'zhujuan'
# last_modify : 2018.02.26

from biocluster.agent import Agent
from biocluster.tool import Tool
from biocluster.core.exceptions import OptionError
from mbio.packages.align.blast.xml2table import xml2table
import pandas as pd
import subprocess, os,re


class AnnoKeggAgent(Agent):
    """
    调用kegg注释需要的脚本 kegg_xlm_mongo.py 获取注释信息以及对ko_name进行取名规则处理
    """

    def __init__(self, parent):
        super(AnnoKeggAgent, self).__init__(parent)
        options = [
            {"name": "kegg_xml", "type": "infile", "format": "align.blast.blast_xml"},  # 输入的比对结果xml文件
            {"name": "anno_kegg", "type": "outfile", "format": "sequence.profile_table"},  # 基因具体的注释信息
            {"name": "level_stat", "type": "outfile", "format": "sequence.profile_table"},  # levvl层级信息
        ]
        self.add_option(options)

    def check_options(self):
        if not self.option("kegg_xml").is_set:
            raise OptionError("必须设置输入文件", code="31202101")
        return True

    def set_resource(self):
        self._cpu = 2
        self._memory = '10G'

    def end(self):
        super(AnnoKeggAgent, self).end()


class AnnoKeggTool(Tool):
    def __init__(self, config):
        super(AnnoKeggTool, self).__init__(config)
        self._version = "1.0"
        self.python_path = self.config.SOFTWARE_DIR + "/miniconda2/bin/python"
        self.python_script = self.config.PACKAGE_DIR + '/bac_comp_genome/kegg_anno.py'

    def run(self):
        """
        运行
        :return:
        """
        super(AnnoKeggTool, self).run()
        self.run_kegg_anno()
        self.set_output()
        self.end()

    def run_kegg_anno(self):
        table = xml2table(self.option('kegg_xml').prop['path'], self.output_dir + '/tmp_kegg_table.xls')
        cmd = '{} {} -i {} -o {}'.format(self.python_path, self.python_script, table, self.output_dir)
        self.logger.info(cmd)
        try:
            subprocess.check_output(cmd, shell=True)
            self.logger.info('运行kegg_anno完成')
        except subprocess.CalledProcessError:
            self.set_error('运行kegg_anno出错')

    def set_output(self):
        pass
