# -*- coding: utf-8 -*-
# __author__ = 'shaohua.yuan'

from biocluster.agent import Agent
from biocluster.tool import Tool
from biocluster.core.exceptions import OptionError
import subprocess
import os
import pandas as pd
import unittest


class CorrelationFdAgent(Agent):
    """
    宏基因组贡献度分析
    author: shaohua.yuan
    last_modify:
    """

    def __init__(self, parent):
        super(CorrelationFdAgent, self).__init__(parent)
        options = [
            {"name": "table1", "type": "infile", "format": "sequence.profile_table", "required": True},
            {"name": "table2", "type": "infile", "format": "sequence.profile_table", "required": True},
            {"name": "coefficient", "type": "string", "default": "spearman"},  # spearman,pearson,kendall
            {"name": "coefficient_value", "type": "float", "default": 0.6},
            {"name": "p_value", "type": "float", "default": 0.05},
            {"name": "trans_t2", "type": "bool", "default": False},  # 表格2是否需要转置
            {"name": "correlation_file", "type": "outfile", "format": "sequence.profile_table"}
        ]
        self.add_option(options)

    def check_options(self):
        if not self.option("coefficient") in ["spearman", "pearson", "kendall"]:
            raise OptionError("相关性系数类型必须为spearman,pearson或kendall", code="")
        return True

    def set_resource(self):
        self._cpu = 3
        self._memory = '15G'

    def end(self):
        super(CorrelationFdAgent, self).end()


class CorrelationFdTool(Tool):
    def __init__(self, config):
        super(CorrelationFdTool, self).__init__(config)
        self.python_path = "/program/Python/bin/python"
        self.script = self.config.PACKAGE_DIR + '/metagenomic/scripts/cal_correlation_fd.py'

    def run(self):
        """
        运行
        :return:
        """
        super(CorrelationFdTool, self).run()
        self.run_network()
        self.set_output()
        self.end()

    def run_network(self):
        profile_table1 = self.option("table1").prop["path"]
        profile_table2 = self.option("table2").prop["path"]
        '''
        table = pd.read_table(profile_table, sep='\t', header=0)
        if len(table) < 2:
            self.set_error("该筛选方法TOP丰度表物种或功能个数小于2无法计算相关性！", code="")
        '''
        cmd = '{} {} -i1 {} -i2 {} -c {} -c_cut {} -p_cut {} -o {}'.format(self.python_path, self.script, profile_table1,
                                                                           profile_table2, self.option('coefficient'),
                                                                           self.option('coefficient_value'),
                                                                           self.option('p_value'),
                                                                           self.output_dir + '/correlation.xls')
        if self.option("trans_t2"):
            cmd += " --trans "
        self.logger.info('开始计算相关性!')
        self.logger.info(cmd)
        command = self.add_command("cal_correlationfd", cmd).run()
        self.wait(command)
        if command.return_code == 0:
            self.logger.info("cal_correlationfd succeed")
        else:
            self.set_error("cal_correlationfd failed", code="")
            raise Exception("cal_correlationfd failed")

    def set_output(self):
        self.option('correlation_file', os.path.join(self.output_dir, 'correlation.xls'))


class TestFunction(unittest.TestCase):
    """
    This is test for the tool. Just run this script to do test.
    """
    def test(self):
        import random
        from mbio.workflows.single import SingleWorkflow
        from biocluster.wsheet import Sheet
        data = {
            "id": "CorrelationFd",
            "type": "tool",
            "name": "meta.association_model.correlation_fd",
            "instant": True,
            "options": dict(
                table1="/mnt/ilustre/users/sanger-dev/workspace/20180925/BetaDiversity_tsg_31835_64490_966855/new_abund_table.xls",
                table2="/mnt/ilustre/users/sanger-dev/workspace/20180925/BetaDiversity_tsg_31835_64490_966855/env_file_input_env.xls",
                trans_t2=True,
            )
        }
        wsheet = Sheet(data=data)
        wf = SingleWorkflow(wsheet)
        wf.run()


if __name__ == '__main__':
    unittest.main()
