# -*- coding: utf-8 -*-
# __author__ = 'shicaiping'
# modified 20210203

from biocluster.agent import Agent
from biocluster.tool import Tool
from biocluster.core.exceptions import OptionError
import os
import unittest
import subprocess
import datetime

starttime = datetime.datetime.now()


class UniqMappedAgent(Agent):
    """
    基因比对接口
    """

    def __init__(self, parent):
        super(UniqMappedAgent, self).__init__(parent)
        options = [
            {"name": "bam_path", "type": "infile", "format": "ref_rna_v2.common"},
        ]
        self.add_option(options)

    def check_options(self):
        if not self.option("bam_path"):
            raise OptionError("请设置bam路径")

    def set_resource(self):
        self._cpu = 2
        self._memory = "10G"

    def end(self):
        super(UniqMappedAgent, self).end()


class UniqMappedTool(Tool):
    def __init__(self, config):
        super(UniqMappedTool, self).__init__(config)
        self.samtools = 'miniconda2/bin/samtools'

    def run_uniq_mapped(self):
        bam_file = self.option("bam_path").prop['path']
        sample = os.path.basename(bam_file).replace(".bam", "")
        cmd1 = "{} view -h {} -O SAM -o {}".format(self.samtools, bam_file, sample + ".sam")
        command1 = self.add_command("bam2sam", cmd1).run()
        self.wait()
        if command1.return_code == 0:
            self.logger.info("bam2sam运行完成")
        else:
            self.set_error("bam2sam运行失败")
        cmd2 = "grep -E \"NH:i:1|^@\" {} > {}".format(sample + ".sam", "uniq.sam")
        try:
            subprocess.check_output(cmd2, shell=True)
            self.logger.info("提取uniq mapped序列完成")
        except:
            self.set_error("提取uniq mapped序列失败")
        out_bam = self.output_dir + "/" + sample + ".bam"
        cmd3 = "{} view -bh {} -o {}".format(self.samtools, "uniq.sam", out_bam)
        command2 = self.add_command("sam2bam", cmd3).run()
        self.wait()
        if command2.return_code == 0:
            self.logger.info("sam2bam运行完成")
        else:
            self.set_error("sam2bam运行失败")

    def run(self):
        super(UniqMappedTool, self).run()
        self.run_uniq_mapped()
        self.end()


class TestFunction(unittest.TestCase):
    def test(self):
        import random
        from mbio.workflows.single import SingleWorkflow
        from biocluster.wsheet import Sheet
        data = {
            'id': 'uniq_mapped_{}_{}'.format(random.randint(1000, 9999), random.randint(1000, 9999)),
            'type': 'tool',
            'name': 'ref_rna_v2.uniq_mapped',
            'instant': False,
            'options': dict(
                bam_path="/mnt/lustre/users/sanger/workspace/20210202/MJ20200603141/RnaseqMapping/Hisat1/output/KF19Q16R_1.bam",
            )
        }
        wsheet = Sheet(data=data)
        wf = SingleWorkflow(wsheet)
        wf.run()

if __name__ == '__main__':
    suite = unittest.TestSuite()
    suite.addTests([TestFunction('test')])
    unittest.TextTestRunner(verbosity=2).run(suite)