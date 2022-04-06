# -*- coding: utf-8 -*-
# __author__ = 'zhouxuan'

from biocluster.agent import Agent
from biocluster.tool import Tool
from biocluster.core.exceptions import OptionError
from mbio.packages.align.blast.xml2table import xml2table
import subprocess
import os


class MgCogAnnoAgent(Agent):
    """
    宏基因cog注释重写，根据比对结果文件从参考库中获取注释信息
    author: zhouxuan
    last_modify: 2017.09.11
    last_modify by: shaohua.yuan
    """
    def __init__(self, parent):
        super(MgCogAnnoAgent, self).__init__(parent)
        options = [
            {"name": "cog_xml", "type": "infile", "format": "align.blast.blast_xml"},  # 比对到string库的xml文件
        ]
        self.add_option(options)
        self._memory_increase_step = 20  # 每次重运行增加内存40G by qingchen.zhang @ 20190527

    def check_options(self):
        if not self.option("cog_xml").is_set:
            raise OptionError("必须设置输入文件", code="31202601")
        return True

    def set_resource(self):
        self._cpu = 2
        self._memory = '10G' ##fix by qingchen.zhang 由5G改为10G

    def end(self):
        result_dir = self.add_upload_dir(self.output_dir)
        result_dir.add_relpath_rules([
            [".", "", "结果输出目录"],
        ])
        super(MgCogAnnoAgent, self).end()


class MgCogAnnoTool(Tool):
    def __init__(self, config):
        super(MgCogAnnoTool, self).__init__(config)
        self._version = "1.0"
        self.python_path = self.config.SOFTWARE_DIR + "/miniconda2/bin/python"
        self.python_script = self.config.PACKAGE_DIR + '/annotation/mg_annotation/mg_cog_annotation.py'

    def run(self):
        """
        运行
        :return:
        """
        super(MgCogAnnoTool, self).run()
        self.run_cog_anno()
        self.set_output()
        self.end()

    def run_cog_anno(self):
        table = xml2table(self.option('cog_xml').prop['path'], self.output_dir + '/tmp_cog_table.xls')
        cmd = '{} {} -i {} -o {}'.format(self.python_path, self.python_script, table,
                                   self.output_dir + "/gene_cog_anno.xls")
        self.logger.info(cmd)
        try:
            subprocess.check_output(cmd, shell=True)
            self.logger.info('运行cog_anno完成')
        except subprocess.CalledProcessError:
            self.set_error('运行cog_anno出错', code="31202601")

    def set_output(self):
        if len(os.listdir(self.output_dir)) == 2:
            self.logger.info("output right")
        else:
            self.logger.info("output wrong")
