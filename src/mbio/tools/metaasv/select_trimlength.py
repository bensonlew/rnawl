# -*- coding: utf-8 -*-
# __author__ = 'qingchen.zhang'

import os
from biocluster.core.exceptions import OptionError
from biocluster.agent import Agent
from biocluster.tool import Tool


class SelectTrimlengthAgent(Agent):
    """
    metaasv 根据输入的fastq文件夹 筛选出最佳的序列长度
    """
    def __init__(self, parent):
        super(SelectTrimlengthAgent, self).__init__(parent)
        options = [
            {"name": "fastq_dir", "type": "infile", "format": "sequence.fastq_dir"},  # 质控序列文件夹
            {"name": "coverage", "type": "int", "default": 90}  # 根据90%的水平进行筛选， 理论上选择最小的文件
        ]
        self.add_option(options)

    def check_options(self):
        """
        检查参数
        :return:
        """
        if not self.option('fastq_dir'):
            raise OptionError('必须输入fastq_dir文件夹')

    def set_resource(self):
        """
        :return:
        """
        self._cpu = 1
        self._memory = "20G"

    def end(self):
        super(SelectTrimlengthAgent, self).end()

class SelectTrimlengthTool(Tool):
    def __init__(self, config):
        super(SelectTrimlengthTool, self).__init__(config)
        self.dir = self.option('fastq_dir').prop['path']
        self.perl_path = "/program/perl-5.24.0/bin/perl"
        self.perl_script = self.config.PACKAGE_DIR + "/metaasv/"
        self.out = self.work_dir + "/all_result.xls"

    def run(self):
        """
        运行
        :return:
        """
        super(SelectTrimlengthTool, self).run()
        self.select_coverage()
        self.set_output()
        self.end()

    def select_coverage(self):
        """
        根据coverage得到每个样本的最短长度
        :return:
        """
        cmd = '{} {}select_coverage.pl {} {} {}'.format(self.perl_path, self.perl_script, self.dir,self.out, self.option("coverage"))
        command = self.add_command('select', cmd).run()
        self.wait(command)
        if command.return_code == 0:
            self.logger.info("select运行完成" )
        else:
            self.set_error("select运行出错!")

    def set_output(self):
        """
        将结果文件复制到output文件夹下面
        :return:
        """
        self.logger.info("设置结果目录")
        if os.path.exists(self.out):
            os.link(self.out, self.output_dir + "/all_result.xls")

