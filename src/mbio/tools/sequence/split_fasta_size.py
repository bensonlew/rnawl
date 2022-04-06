# -*- coding: utf-8 -*-
# __author__ = 'zouguanqing'
from biocluster.agent import Agent
from biocluster.tool import Tool
from biocluster.core.exceptions import OptionError



class SplitFastaSizeAgent(Agent):
    """
    SplitFastaSize:将fasta文件按文件大小拆分
    version 1.0
    author: zouguanqing
    """

    def __init__(self, parent):
        super(SplitFastaSizeAgent, self).__init__(parent)
        options = [
            {"name": "fasta", "type": "infile", "format": "sequence.fasta"},
            {"name": "unit_size", "type": "int", "default": 100000000},  # 序列数
        ]
        self.add_option(options)


    def check_options(self):
        """
        检查参数是否正确
        """
        if not self.option("fasta").is_set:
            raise OptionError("请传入fasta序列文件0", code="34003201")
        if not isinstance(self.option('unit_size'), int):
            raise OptionError("文件大小必须为整数", code="34003202")


    def set_resource(self):
        """
        所需资源
        """
        self._cpu = 2
        self._memory = '20G'


class SplitFastaSizeTool(Tool):
    """
    version 1.0
    """

    def __init__(self, config):
        super(SplitFastaSizeTool, self).__init__(config)
        self.python = 'miniconda2/bin/python'
        self.divide_src = self.config.PACKAGE_DIR + '/sequence/fasta_split_by_size.py'

    def split_fasta(self):
        """
        """
        fasta_file = self.option('fasta').prop['path']
        unit_size = self.option('unit_size')
        cmd = '{} {} {} {} size'.format(self.python, self.divide_src, fasta_file, unit_size)
        self.logger.info(cmd)
        command = self.add_command('divide', cmd)
        command.run()
        if command.return_code == 0:
            self.logger.info("运行divide完成")
        else:
            self.set_error("运行出错!", code="34003201")

    def run(self):
        super(SplitFastaSizeTool, self).run()
        self.split_fasta()
        self.end()
