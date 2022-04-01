# -*- coding: utf-8 -*-
# __author__ = 'shaohua.yuan'

from biocluster.agent import Agent
from biocluster.tool import Tool
from biocluster.core.exceptions import OptionError
import subprocess
import os
import pandas as pd


class NetworkCorAgent(Agent):
    """
    宏基因组相关性分析
    author: shaohua.yuan
    last_modify:
    """

    def __init__(self, parent):
        super(NetworkCorAgent, self).__init__(parent)
        options = [
            {"name": "correlation_file", "type": "infile", "format": "sequence.profile_table"},
        ]
        self.add_option(options)

    def check_options(self):
        self.logger.info(self.option("correlation_file").prop["path"])
        if not self.option("correlation_file").is_set:
            raise OptionError("必须设置输入相关性文件", code="32701401")
        return True

    def set_resource(self):
        self._cpu = 2
        self._memory = '15G'

    def end(self):
        super(NetworkCorAgent, self).end()


class NetworkCorTool(Tool):
    def __init__(self, config):
        super(NetworkCorTool, self).__init__(config)
        self.python_path = "/program/Python/bin/python"
        self.cmd_path = self.config.PACKAGE_DIR + '/metagenomic/scripts/corr_net_calc.py'

    def run(self):
        """
        运行
        :return:
        """
        super(NetworkCorTool, self).run()
        self.run_network()
        self.set_output()
        self.end()

    def run_network(self):
        correlation_file = self.option("correlation_file").prop["path"]
        table = pd.read_table(correlation_file , sep='\t', header=0)
        if len(table) < 1:
            self.set_error("该筛选阈值筛选的物种或功能小于2无法计算网络相关系数！", code="32701401")
        cmd = '{} {} -i {} -o {}'.format(self.python_path, self.cmd_path, correlation_file, self.output_dir)
        self.logger.info('开始生成网络并进行计算')
        command = self.add_command("network_cor", cmd).run()
        self.logger.info(cmd)
        self.wait(command)
        if command.return_code == 0:
            self.logger.info("相关性网络生成成功")
        else:
            self.set_error("相关性网络生成失败", code="32701402")
            raise Exception("相关性网络生成失败")

    def set_output(self):
        self.logger.info("Network output设置")
        if len(os.listdir(self.output_dir)) == 8 :
            self.logger.info("Network output right")
        else:
            self.logger.info("Network output wrong")