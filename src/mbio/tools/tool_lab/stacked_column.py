# -*- coding: utf-8 -*-
# __author__ = 'fwy'
import os
import shutil
from biocluster.core.exceptions import OptionError
from biocluster.agent import Agent
from biocluster.tool import Tool
from biocluster.config import Config
import re
import pandas as pd

import unittest


class StackedColumnAgent(Agent):
    """
    Used for fasta seq stat .
    """
    def __init__(self, parent):
        super(StackedColumnAgent, self).__init__(parent)
        options = [
            {"name": "input_file", "type": "infile", "format": "ref_rna_v2.common"},  # 输入文件，可以是表达量等文件
            {"name": "filter_value", "type": "float", "default":0.05},  # 过滤阈值，低于该阈值的归为others
            {"name": "sep", "type": "string", "default":"tab"},  #表格分隔符
        ]
        self.add_option(options)
        self.step.add_steps("stacked_column")
        self.on('start', self.stepstart)
        self.on('end', self.stepfinish)

    def stepstart(self):
        self.step.stacked_column.start()
        self.step.update()

    def stepfinish(self):
        self.step.stacked_column.finish()
        self.step.update()

    def check_options(self):
        """
        重写参数检测函数
        :return:
        """
        if not self.option('input_file').is_set:
            raise OptionError('必须输入待统计表格')
        return True

    def set_resource(self):
        """
        设置所需资源，需在之类中重写此方法 self._cpu ,self._memory
        :return:
        """
        self._cpu = 2
        self._memory = "5G"

    def end(self):
        result_dir = self.add_upload_dir(self.output_dir)
        super(StackedColumnAgent, self).end()


class StackedColumnTool(Tool):
    def __init__(self, config):
        super(StackedColumnTool, self).__init__(config)
        self._version = "v1.0.1"
        self.python_path = '/program/Python/bin/python'
        self.perl = 'program/perl/perls/perl-5.24.0/bin/perl'
        self.tool_path=self.config.PACKAGE_DIR+"/tool_lab/stacked_column.py"

    def run(self):
        """
        运行
        :return:
        """
        super(StackedColumnTool, self).run()
        self.run_stacked_column()
        self.set_output()
        self.end()

    def run_stacked_column(self):
        self.logger.info("开始对fasta文件进行统计")
        stat_cmd = '{} {} '.format(self.python_path,self.tool_path)
        stat_cmd += '-i {} '.format(self.option("input_file").prop["path"])
        stat_cmd += '-s {} '.format(self.option("sep"))
        stat_cmd += '-v {} '.format(self.option("filter_value"))
        stat_cmd += '-o stacked_column.xls '
        command = self.add_command("stacked_column", stat_cmd, ignore_error=True)
        command.run()
        self.wait(command)
        if command.return_code == 0:
            self.logger.info("运行stacked_column完成")
        else:
            self.set_error("运行stacked_column运行出错!")
            return False



    def set_output(self):
        """
        将结果文件复制到output文件夹下面
        :return:
        """
        self.logger.info("设置结果目录")
        try:
                result_detail=os.path.join(self.work_dir, "stacked_column.xls")
                link_detail = os.path.join(self.output_dir, "stacked_column.xls")
                if os.path.exists(link_detail):
                    os.remove(link_detail)
                os.link(result_detail, link_detail)
        except Exception as e:
            self.logger.info("设置结果目录失败{}".format(e))


class TestFunction(unittest.TestCase):
    """
    This is test for the tool. Just run this script to do test.
    """
    def test(self):
        import random
        from mbio.workflows.single import SingleWorkflow
        from biocluster.wsheet import Sheet
        test_dir = "/mnt/ilustre/users/sanger-dev/app/database/Genome_DB_finish/fungi/Saccharomyces_cerevisiae/Ensemble_release_39/"
        data = {
            "id": "seq_stat" + str(random.randint(1, 10000)),
            "type": "tool",
            "name": "tool_lab.stacked_column",
            "instant": True,
            "options": dict(
                # fasta_file=test_dir + "/" + "dna/Saccharomyces_cerevisiae.dna.toplevel.fa",
                input_file ="/mnt/ilustre/users/sanger-dev/sg-users/fuwenyao/tool_lab/stacked_column/in_otu_table.xls" ,
                filter_value=0.3,
                sep="tab"
                # step=500
            )
        }
        wsheet = Sheet(data=data)
        wf = SingleWorkflow(wsheet)
        wf.run()


if __name__ == '__main__':
    unittest.main()
