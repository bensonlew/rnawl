# -*- coding: utf-8 -*-
# __author__ = 'zhangyitong'

from biocluster.agent import Agent
from biocluster.tool import Tool
from biocluster.core.exceptions import OptionError
import os
import unittest
import re


class ReadsRepairAgent(Agent):
    """
    群体进化，vcftool过滤vcf文件
    """
    def __init__(self, parent):
        super(ReadsRepairAgent, self).__init__(parent)
        options = [
            {"name": "r1", "type": "infile", "format": "ref_rna_v2.fastq"},
            {'name': 'r2', "type": 'infile', 'format': 'ref_rna_v2.fastq'},
        ]
        self.add_option(options)
        self.step.add_steps('reads_repair')
        self.on('start', self.step_start)
        self.on('end', self.step_end)

    def step_start(self):
        self.step.reads_repair.start()
        self.step.update()

    def step_end(self):
        self.step.reads_repair.finish()
        self.step.update()

    def check_options(self):
        if not self.option("r1").is_set:
            raise OptionError("请上传R1序列")
        if not self.option("r2").is_set:
            raise OptionError("请上传R2序列")
        return True

    def set_resource(self):
        """
        """
        self._cpu = 1
        self._memory = "30G"

    def end(self):
        super(ReadsRepairAgent, self).end()


class ReadsRepairTool(Tool):
    def __init__(self, config):
        super(ReadsRepairTool, self).__init__(config)
        java_dir = os.path.join(self.config.SOFTWARE_DIR, 'program/sun_jdk1.8.0/bin')
        self.set_environ(PATH=java_dir)
        self.program = {
            'shell': 'program/sh',
        }
        self.script = {
            'repair':  os.path.join(self.config.SOFTWARE_DIR, 'bioinfo/tool_lab/bbmap/repair.sh'),
            'reformat': os.path.join(self.config.SOFTWARE_DIR, 'bioinfo/tool_lab/bbmap/reformat.sh'),
        }
        self.file = {
            'out': os.path.join(self.work_dir, 'mid.fq.gz'),
            'outsingle': os.path.join(self.work_dir, 'remove_reads.fq.gz')
        }

    def run_repair(self):
        """
        """
        cmd = "{} {} ".format(self.program['shell'], self.script['repair'])
        cmd += "in1={} ".format(self.option('r1').prop['path'])
        cmd += 'in2={} '.format(self.option('r2').prop['path'])
        cmd += "out={} ".format(self.file['out'])
        cmd += 'outsingle={} overwrite=t'.format(self.file['outsingle'])
        self.logger.info(cmd)
        self.logger.info("开始进行bbmap repair")
        command = self.add_command("bbmap_repair", cmd).run()  # 必须小写，
        self.wait(command)
        if command.return_code == 0:
            self.logger.info("bbmap repair完成！")
        else:
            self.set_error("bbmap repair出错！")
        self.run_reformat()

    def run_reformat(self):
        out1 = os.path.join(self.output_dir, 'repair_' + os.path.basename(self.option('r1').prop['path']))
        out2 = os.path.join(self.output_dir, 'repair_' + os.path.basename(self.option('r2').prop['path']))
        cmd = "{} {} ".format(self.program['shell'], self.script['reformat'])
        cmd += "in={} ".format(self.file['out'])
        cmd += 'out1={} '.format(out1)
        cmd += "out2={} ".format(out2)
        cmd += 'addcolon overwrite=t'
        self.logger.info(cmd)
        self.logger.info("开始进行bbmap reformat")
        command = self.add_command("bbmap_reformat", cmd).run()  # 必须小写，
        self.wait(command)
        if command.return_code == 0:
            self.logger.info("bbmap reformat运行完成！")
        else:
            self.set_error("bbmap reformat运行出错！")

    def run(self):
        super(ReadsRepairTool, self).run()
        self.run_repair()
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