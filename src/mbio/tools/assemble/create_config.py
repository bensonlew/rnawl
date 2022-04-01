# -*- coding: utf-8 -*-
# __author__ = 'hao.gao'

import os, sys
# import shutil
from biocluster.core.exceptions import OptionError
from biocluster.agent import Agent
from biocluster.tool import Tool


class CreateConfigAgent(Agent):
    """
    细菌基因组SOAPdenovo2组装config文件配置
    version: v1
    author: hao.gao
    last_modify: 2018.01.25
    """

    def __init__(self, parent):
        super(CreateConfigAgent, self).__init__(parent)
        options = [
            {"name": "PE_list", "type": "infile", "format": "meta.otu.otu_table"},  # 样品的PE文库信息表
            {"name": "MP_list", "type": "infile", "format": "meta.otu.otu_table"},  # 样品的MP文库信息表
            {"name": "config_file", "type": "outfile", "format": "bacgenome.config_file"},   #配置好的config文件
            {'name': 'sample_name', "type": "string"},  # 样本名
        ]
        self.add_option(options)


    def check_options(self):
        """
        检查参数
        :return:
        """
        if not self.option('PE_list').is_set:
            raise OptionError('必须输入PE_list文件', code="31300301")
        if not self.option('sample_name'):
            raise OptionError('必须输入sample_name文件', code="31300302")

    def set_resource(self):
        """
        :return:
        """
        self._cpu = 2
        self._memory = "5G"

    def end(self):
        super(CreateConfigAgent, self).end()


class CreateConfigTool(Tool):
    def __init__(self, config):
        super(CreateConfigTool, self).__init__(config)
        self.pe_list = self.option('PE_list').prop['path']
        self.sample_name = self.option('sample_name')
        self.python_path = "/program/Python/bin/python"
        self.python_script = self.config.PACKAGE_DIR + "/bacgenome/"
        self.config_file = self.work_dir + "/" + self.sample_name + ".config"

    def run(self):
        """
        运行
        :return:
        """
        super(CreateConfigTool, self).run()
        self.init_config()
        self.set_output()
        self.end()

    def init_config(self):
        if self.option('MP_list').is_set:
            self.mp_list = self.option('MP_list').prop['path']
            cmd = '{} {}config_file.py -p {} -m {} -o {}'.format(self.python_path, self.python_script, self.pe_list, self.mp_list, self.config_file)
            command = self.add_command('pe_mp_config', cmd).run()
            self.wait(command)
            if command.return_code == 0:
                self.logger.info("pe_mp_config运行完成" )
            else:
                self.set_error("pe_mp_config运行出错!", code="31300301")
        else:
            cmd = '{} {}config_file.py -p {} -o {}'.format(self.python_path, self.python_script, self.pe_list, self.config_file)
            command = self.add_command('pe_config', cmd).run()
            self.wait(command)
            if command.return_code == 0:
                self.logger.info("pe_config运行完成")
            else:
                self.set_error("pe_config运行出错!", code="31300302")


    def set_output(self):
        """
        将结果文件复制到output文件夹下面
        :return:
        """
        self.logger.info("设置结果目录")
        self.option('config_file').set_path(self.config_file)
        self.logger.info("设置SOAPdenovo2分析结果目录成功")
