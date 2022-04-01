# -*- coding: utf-8 -*-


import os, shutil
from biocluster.core.exceptions import OptionError
from biocluster.agent import Agent
from biocluster.tool import Tool
import unittest



class VarianceHomogAgent(Agent):


    def __init__(self, parent):
        super(VarianceHomogAgent, self).__init__(parent)
        options = [
            {"name": "table", "type": "infile", "format": "tool_lab.simple"},
            {"name": "group", "type": "infile", "format": "tool_lab.simple"},
            {"name": "method", "type": "string", "default": "bartlett"}, #bartlett\levene\ligner-killeen
        ]
        self.add_option(options)

    def check_options(self):
        """
        检查参数
        :return:
        """
        if not self.option('table').is_set or  not self.option('group').is_set:
            raise OptionError('must input table and group file')

        if self.option('method') not in ['bartlett','levene','ligner-killeen']:
            raise OptionError("method must in ['bartlett','levene','ligner-killeen']")

        return True

    def set_resource(self):
        """
        设置所需资源，需在子类中重写此方法 self._cpu ,self._memory
        :return:
        """
        self._cpu = 1
        self._memory = "5G"

    def end(self):
        super(VarianceHomogAgent, self).end()


class VarianceHomogTool(Tool):
    def __init__(self, config):
        super(VarianceHomogTool, self).__init__(config)
        self.python_path =  "/program/Python/bin/python"
        self.script = self.config.PACKAGE_DIR + '/tool_lab/variance_homog.py'

    def run(self):
        """
        运行
        :return:
        """
        super(VarianceHomogTool, self).run()
        self.logger.info("开始运行命令！")
        self.start_cmd_fun()

        self.set_output()

    def start_cmd_fun(self):
        if self.option('method') == 'ligner-killeen':
            outfile = self.output_dir + '/Homogeneity_Variances_f'+self.option('method')+'.xls'
        else:
            outfile = self.output_dir + '/Homogeneity_Variances_'+self.option('method')+'.xls'
        cmd = self.python_path + ' {} {}  {} {} {}'.format(
            self.script, self.option('table').path, self.option('group').path,self.option('method'), outfile)

        self.logger.info(cmd)
        command = self.add_command('variance_homog', cmd, ignore_error=True).run()
        self.wait(command)
        if command.return_code == 0:
            self.logger.info("run succeed")
        else:
            self.set_error("run failed")
            raise Exception("run failed")

    def set_output(self):
        """
        将结果文件复制到output文件夹下面
        :return:
        """

        self.end()


