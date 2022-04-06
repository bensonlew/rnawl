# -*- coding: utf-8 -*-
# __author__ = 'zouxuan'
from biocluster.agent import Agent
from biocluster.tool import Tool
from biocluster.core.exceptions import OptionError
import os
import re


class AddGeneInfoAgent(Agent):
    """
    
    author: zouxuan
    last_modify: 20180403
    """

    def __init__(self, parent):
        super(AddGeneInfoAgent, self).__init__(parent)
        options = [
            {"name": "sequence", "type": "infile", "format": "sequence.fasta"},  # 序列文件
            {"name": "sample", "type": "string"},  # 样品名
            {"name": "table", "type": "infile", "format": "sequence.profile_table"},  # 待修改的二维表格
            {"name": "split", "type": "bool","default": False}  # 是否按location对文件切割
        ]
        self.add_option(options)

    def check_options(self):
        if not self.option("sequence").is_set:
            raise OptionError("必须设置输入序列文件", code="31400101")
        if not self.option("table").is_set:
            raise OptionError("必须设置待修改的二维表格", code="31400102")
        return True

    def set_resource(self):
        self._cpu = 1
        self._memory = '3G'

    def end(self):
        result_dir = self.add_upload_dir(self.output_dir)
        result_dir.add_relpath_rules([
            [".", "", "结果输出目录"],
        ])
        super(AddGeneInfoAgent, self).end()


class AddGeneInfoTool(Tool):
    def __init__(self, config):
        super(AddGeneInfoTool, self).__init__(config)
        self._version = "1.0"
        self.script_path = self.config.PACKAGE_DIR + '/bacgenome/add_gene_info.py '
        self.python_path = '/miniconda2/bin/python'

    def run(self):
        """
        运行
        :return:
        """
        super(AddGeneInfoTool, self).run()
        self.get_file()
        self.end()

    def get_file(self):
        cmd = '{} {} -i {} -f {} -s {} -d {}'.format(self.python_path, self.script_path,
                                         self.option('table').prop['path'],self.option('sequence').prop['path'],
                                                self.option('sample'),self.option('split'))
        command = self.add_command('get_file', cmd).run()
        self.wait(command)
        if command.return_code == 0:
            self.logger.info("文件生成成功")
        else:
            self.set_error("文件生成失败", code="31400101")
            self.set_error("文件生成失败", code="31400102")

