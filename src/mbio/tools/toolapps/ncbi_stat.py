# -*- coding: utf-8 -*-
# __author__ = 'hao.gao'

from biocluster.core.exceptions import OptionError
from biocluster.agent import Agent
from biocluster.tool import Tool
import os
import unittest


class NcbiStatAgent(Agent):
    '''
    NCBI下载统计16s和看家基因的数据
    '''
    def __init__(self, parent):
        super(NcbiStatAgent, self).__init__(parent)
        options = [
            {"name": "sample_list", "type": "string"},  # 样品名称，sample1;sample2;...
            {'name': 'house_dir', 'type': 'string'},  ## 看家基因的结果目录
            {'name': 's16', 'type': 'string'},  ## 16结果统计文件
        ]
        self.add_option(options)

    def check_options(self):
        """
        参数二次检查
        :return:
        """
        if not self.option('sample_list'):
            raise OptionError("必须输入sample_list！")
        return True

    def set_resource(self):
        """
        设置资源
        :return:
        """
        self._cpu = '2'
        self._memory = '5G'

    def end(self):
        """
        运行结束
        :return:
        """
        super(NcbiStatAgent, self).end()

class NcbiStatTool(Tool):
    """
    NCBI下载统计16s和看家基因的数据
    """
    def __init__(self, config):
        super(NcbiStatTool, self).__init__(config)
        self.python = "/miniconda2/bin/python"
        self.script = self.config.PACKAGE_DIR + "/toolapps/"


    def run(self):
        """
        运行
        :return:
        """
        self.logger.info("开始运行tool!")
        super(NcbiStatTool, self).run()
        self.run_stat()
        self.set_output()
        self.end()

    def run_stat(self):
        """
        运行gtdbtk的sh，里面包含分析的三步
        :return:
        """
        self.logger.info("分析开始!")
        cmd = '{} {}ncbi_stat.py -i \"{}\" -s {} -d {} -o {}'.format(self.python, self.script, self.option("sample_list"), self.option("s16"), self.option("house_dir"), self.output_dir + "/all_database.stat.xls")
        self.logger.info(cmd)
        command = self.add_command("run_stat", cmd)
        command.run()
        self.wait()
        if command.return_code == 0:
            self.logger.info('run_stat软件运行成功！')
        else:
            self.set_error('run_stat软件运行失败！')

    def set_output(self):
        """
        设置结果文件目录
        :return:
        """
        self.logger.info('开始设置结果文件目录')