# -*- coding: utf-8 -*-
# __author__ = 'shaohua.yuan'

from biocluster.agent import Agent
from biocluster.tool import Tool
from biocluster.core.exceptions import OptionError
import subprocess
import os
import pandas as pd

class CorrelationAgent(Agent):
    """
    宏基因组贡献度分析
    author: shaohua.yuan
    last_modify:
    """

    def __init__(self, parent):
        super(CorrelationAgent, self).__init__(parent)
        options = [
            {"name": "profile_table", "type": "infile", "format": "sequence.profile_table"},
            {"name": "coefficient", "type": "string", "format": "spearman"},  # spearman,pearson,kendall
            {"name": "coefficient_value", "type": "float", "default": 0.6},
            {"name": "p_value", "type": "float", "default": 0.05},
            {"name": "correlation_file", "type": "outfile", "format": "sequence.profile_table"}
        ]
        self.add_option(options)

    def check_options(self):
        if not self.option("profile_table").is_set:
            raise OptionError("必须设置输入丰度文件", code="32701101")
        return True

    def set_resource(self):
        self._cpu = 3
        self._memory = '15G'

    def end(self):
        super(CorrelationAgent, self).end()


class CorrelationTool(Tool):
    def __init__(self, config):
        super(CorrelationTool, self).__init__(config)
        self.python_path = "/miniconda2/bin/python"
        self.script = self.config.PACKAGE_DIR + '/metagenomic/scripts/cal_correlation.py'

    def run(self):
        """
        运行
        :return:
        """
        super(CorrelationTool, self).run()
        self.run_network()
        self.set_output()
        self.end()

    def run_network(self):
        profile_table = self.option("profile_table").prop["path"]
        table = pd.read_table(profile_table , sep='\t', header=0)
        if len(table) < 2:
            self.set_error("该筛选方法TOP丰度表物种或功能个数小于2无法计算相关性！", code="32701101")
        cmd = '{} {} -i {} -c {} -c_cut {} -p_cut {} -o {}'.format(self.python_path, self.script, profile_table,
                                                                   self.option('coefficient'),
                                                                   self.option('coefficient_value'),
                                                                   self.option('p_value'),
                                                                   self.output_dir + '/correlation.xls')
        self.logger.info('开始计算相关性!')
        self.logger.info(cmd)
        command = self.add_command("cal_correlation", cmd).run()
        self.wait(command)
        if command.return_code == 0:
            self.logger.info("cal_correlation succeed")
        else:
            self.set_error("cal_correlation failed", code="32701102")
            raise Exception("cal_correlation failed")

    def set_output(self):
        self.option('correlation_file', os.path.join(self.output_dir, 'correlation.xls'))
