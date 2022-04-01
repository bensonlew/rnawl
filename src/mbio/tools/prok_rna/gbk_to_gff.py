# -*- coding: utf-8 -*-
# __author__ = 'zhangyitong'
from biocluster.agent import Agent
from biocluster.tool import Tool
from biocluster.core.exceptions import OptionError
import os
import re
import ConfigParser
import unittest
import glob


class GbkToGffAgent(Agent):
    """
    已知miRNA鉴定
    """
    def __init__(self, parent):
        super(GbkToGffAgent, self).__init__(parent)
        options = [
            {"name": "gbk_file", "type": "infile", "format": "prok_rna.common"},  # 输入gbk文件
            {'name': 'gff_file', 'type': 'outfile', 'format': 'prok_rna.common'}
        ]
        self.add_option(options)
        self.step.add_steps("gbk2gff")
        self.on('start', self.stepstart)
        self.on('end', self.stepfinish)

    def stepstart(self):
        self.step.gbk2gff.start()
        self.step.update()

    def stepfinish(self):
        self.step.gbk2gff.finish()
        self.step.update()

    def check_options(self):
        """
        重写参数检测函数
        :return:
        """
        if not self.option("gbk_file").is_set:
            raise OptionError("必须提供输入文件", code = "35002401")
        self.set_resource()
        return True

    def set_resource(self):
        """
        设置所需资源，需在类中重写此方法 self._cpu ,self._memory
        :return:
        """
        self._cpu = 1
        self._memory = '5G'

    def end(self):
        super(GbkToGffAgent, self).end()


class GbkToGffTool(Tool):
    def __init__(self, config):
        super(GbkToGffTool, self).__init__(config)
        self.set_environ(PATH=os.path.join(self.config.SOFTWARE_DIR, 'program/Python/bin'))
        self.program = {
            'python': 'program/Python/bin/python',
        }
        self.script = {
            'gbk2gff': os.path.join(self.config.PACKAGE_DIR, "prok_rna/gbk_to_gff.py"),
        }
        self.file = {
            'gff_file': os.path.join(self.work_dir, os.path.basename(self.option('gbk_file').prop['path']) + '.gff')
        }

    def run(self):
        """
        运行
        :return:
        """
        super(GbkToGffTool, self).run()
        self.run_gbk2gff()
        self.set_output()
        self.end()

    def run_gbk2gff(self):
        self.logger.info('转换gbk格式文件至gff')
        cmd = "{} {} {}".format(self.program['python'], self.script['gbk2gff'], self.option('gbk_file').prop['path'])
        cmd += " {}".format(self.file['gff_file'])
        self.logger.info("开始进行gbk文件转换")
        command = self.add_command("gbk_to_gff", cmd)
        command.run()
        self.wait()
        if command.return_code == 0:
            self.logger.info("gbk文件转换完成!")
        else:
            self.set_error("gbk文件转换出错！")

    def set_output(self):
        self.option("gff_file", self.file['gff_file'])


class TestFunction(unittest.TestCase):

    def test(self):
        import random
        from mbio.workflows.single import SingleWorkflow
        from biocluster.wsheet import Sheet
        data = {
            "id": "InfernalRfam" + str(random.randint(1, 10000)),
            "type": "tool",
            "name": "prok_rna.infernal_rfam",
            "instant": False,
            "options": dict(
                query="/mnt/ilustre/users/sanger-dev/workspace/20210602/Single_Srna_3999/Srna2/Rockhopper__1/Rockhopper_Results/genome.predicted_RNA.fa",
                evalue='1e-5'
                #config="/mnt/ilustre/users/sanger-dev/sg-users/shicaiping/miRNA/quantifier_test/Uniq.cfg.ini"
            )
        }
        wsheet = Sheet(data=data)
        wf = SingleWorkflow(wsheet)
        wf.run()

if __name__ == '__main__':
    unittest.main()
