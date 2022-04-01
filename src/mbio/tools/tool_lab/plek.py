# -*- coding: utf-8 -*-
# __author__ = 'zhangyitong'

from biocluster.agent import Agent
from biocluster.tool import Tool
from biocluster.core.exceptions import OptionError
import os
import unittest
import glob
import re


class PlekAgent(Agent):
    """
    群体进化，vcftool过滤vcf文件
    """
    def __init__(self, parent):
        super(PlekAgent, self).__init__(parent)
        options = [
            {"name": "infasta", "type": "infile", "format": "ref_rna_v2.fasta"},
        ]
        self.add_option(options)
        self.step.add_steps('plek')
        self.on('start', self.step_start)
        self.on('end', self.step_end)

    def step_start(self):
        self.step.plek.start()
        self.step.update()

    def step_end(self):
        self.step.plek.finish()
        self.step.update()

    def check_options(self):
        if not self.option("infasta").is_set:
            raise OptionError("必须设置输入fasta格式的lncRNA序列文件")
        return True

    def set_resource(self):
        """
        """
        self._cpu = 8
        self._memory = "20G"

    def end(self):
        super(PlekAgent, self).end()


class PlekTool(Tool):
    def __init__(self, config):
        super(PlekTool, self).__init__(config)
        self.set_environ(PATH=os.path.join(self.config.SOFTWARE_DIR, 'program/Python/bin'))
        self.program = {
            'python': 'program/Python/bin/python',
        }
        self.script = {
            'plek': os.path.join(self.config.SOFTWARE_DIR, 'bioinfo/tool_lab/PLEK.1.2/PLEK.py'),
        }

    def run_plek(self):
        """
        """
        cmd = "{} {} ".format(self.program['python'], self.script['plek'])
        cmd += "-fasta {} ".format(self.option('infasta').prop['path'])
        cmd += '-out {} '.format('predicted_result')
        cmd += "-thread 8"
        # cmd += '-isoutmsg 1 -isrmtempfile 0'
        self.logger.info(cmd)
        self.logger.info("开始运行plek预测lncRNA")
        command = self.add_command("run_plek", cmd).run()  # 必须小写，
        self.wait(command)
        if command.return_code == 0:
            self.logger.info("plek预测lncRNA完成！")
        else:
            self.set_error("plek预测lncRNA出错！")

    def set_output(self):
        o_file = glob.glob(os.path.join(self.work_dir, 'predicted_result*'))
        if len(o_file) > 0:
            os.link(o_file[0], os.path.join(self.output_dir, os.path.basename(o_file[0])))

    def run(self):
        super(PlekTool, self).run()
        self.run_plek()
        self.set_output()
        self.end()


class TestFunction(unittest.TestCase):
    """
    This is test for the tool. Just run script to do test.
    """

    def test(self):
        import random
        from mbio.workflows.single import SingleWorkflow
        from biocluster.wsheet import Sheet
        data = {
            "id": "gwas_plink_{}_{}".format(random.randint(1000, 9999), random.randint(1000, 9999)),
            "type": "tool",
            "name": "tool_lab.gwas.plink",
            "instant": False,
            "options": dict(
                vcf_file='/mnt/ilustre/users/sanger-dev/workspace/20210519/Single_gwas_filter_3931_6429/VcftoolsFilter/output/pop.recode.vcf',
            )
        }
        wsheet = Sheet(data=data)
        wf = SingleWorkflow(wsheet)
        wf.run()


if __name__ == "__main__":
    suite = unittest.TestSuite()
    suite.addTests([TestFunction("test")])
    unittest.TextTestRunner(verbosity=2).run(suite)