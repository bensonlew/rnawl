# -*- coding: utf-8 -*-
# __author__ = 'zhangyitong'

from biocluster.agent import Agent
from biocluster.tool import Tool
from biocluster.core.exceptions import OptionError
import os
import unittest
import re


class BsaCalcAgent(Agent):
    """
    """
    def __init__(self, parent):
        super(BsaCalcAgent, self).__init__(parent)
        options = [
            {"name": "vcf_file", "type": "infile", "format": "dna_evolution.vcf"},
            {"name": "wp", "type": "string"},   # wild parent
            {"name": "mp", "type": "string"},       # mutant parent
            {"name": "wb", "type": "string"},  # wild bulk
            {"name": "mb", "type": "string"},  # mutant bulk
            {"name": "pop_index", "type": "outfile", "format": "ref_rna_v2.common"},
        ]
        self.add_option(options)
        self.step.add_steps('bsa_calc')
        self.on('start', self.step_start)
        self.on('end', self.step_end)

    def step_start(self):
        self.step.bsa_calc.start()
        self.step.update()

    def step_end(self):
        self.step.bsa_calc.finish()
        self.step.update()

    def check_options(self):
        if not self.option("vcf_file").is_set:
            raise OptionError("请设置vcf_file")
        if not self.option('wb'):
            raise OptionError("请设置输入野生型混池信息")
        if not self.option('mb'):
            raise OptionError("请设置输入突变型混池信息")
        return True

    def set_resource(self):
        """
        运行所需资源
        """
        self._cpu = 1
        self._memory = "20G"

    def end(self):
        super(BsaCalcAgent, self).end()


class BsaCalcTool(Tool):
    def __init__(self, config):
        super(BsaCalcTool, self).__init__(config)
        self.program = {
            'perl': 'miniconda2/bin/perl',
        }
        self.script = {
            'vcf2table': os.path.join(self.config.PACKAGE_DIR, 'tool_lab/bsa/vcf2table.pl'),
            'bsa_calc': os.path.join(self.config.PACKAGE_DIR, 'tool_lab/bsa/bsa_calc.pl')
        }
        self.file = {
            'pop_table': os.path.join(self.output_dir, 'pop.table'),
            'pop_index': os.path.join(self.output_dir, 'pop.index')
        }

    def vcf2table(self):
        """
        """
        cmd = "{} {} ".format(self.program['perl'], self.script['vcf2table'])
        cmd += '-vcf {} '.format(self.option('vcf_file').prop['path'])
        cmd += '-out {} '.format(self.file['pop_table'])
        cmd += '-pldep 0 -phdep 0 -bldep 0 -bhdep 1000 -popt F2 '
        if self.option('wp') and self.option('mp'):
            cmd += '-wp {} -mp {} '.format(self.option('wp'), self.option('mp'))
        cmd += '-wb {} -mb {}'.format(self.option('wb'), self.option('mb'))
        self.logger.info(cmd)
        self.logger.info("开始进行vcf2table")
        command = self.add_command("vcf2table", cmd, ignore_error=True)
        command.run()  # 必须小写，
        self.wait(command)
        if command.return_code == 0:
            self.logger.info("vcf2table完成！")
        else:
            self.set_error("vcf2table出错！")
        with open(self.file['pop_table'], 'r') as pop:
            pop.readline()
            if pop.readline():
                self.bsa_calc()

    def bsa_calc(self):
        """
        """
        cmd = "{} {} ".format(self.program['perl'], self.script['bsa_calc'])
        cmd += '-table {} '.format(self.file['pop_table'])
        cmd += '-out {} '.format(self.file['pop_index'])
        if self.option('wp') and self.option('mp'):
            cmd += '-wp {} -mp {} '.format(self.option('wp'), self.option('mp'))
        cmd += '-wb {} -mb {} '.format(self.option('wb'), self.option('mb'))
        cmd += '-popt F2'
        self.logger.info(cmd)
        self.logger.info("开始进行bsa calculation")
        command = self.add_command("bsa_calc", cmd).run()  # 必须小写，
        self.wait(command)
        if command.return_code == 0:
            self.logger.info("bsa calculation完成！")
        else:
            self.set_error("bsa calculation出错！")
        if os.path.exists(self.file['pop_index']):
            self.logger.info("设置pop_index成功")
            self.option("pop_index", self.file['pop_index'])

    def run(self):
        super(BsaCalcTool, self).run()
        self.vcf2table()
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
            "id": "bsa_bsa_calc_{}_{}".format(random.randint(1000, 9999), random.randint(1000, 9999)),
            "type": "tool",
            "name": "tool_lab.bsa.bsa_calc",
            "instant": False,
            "options": dict(
                vcf_file='/mnt/ilustre/users/sanger-dev/sg-users/zhangyitong/test/tool_052021/bsa/pop.final.vcf',
                # wp='',
                # mp='HQS1',
                mb='HQS1',
                wb='F44_mix',
            )
        }
        wsheet = Sheet(data=data)
        wf = SingleWorkflow(wsheet)
        wf.run()


if __name__ == "__main__":
    suite = unittest.TestSuite()
    suite.addTests([TestFunction("test")])
    unittest.TextTestRunner(verbosity=2).run(suite)