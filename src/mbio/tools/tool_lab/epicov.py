# -*- coding: utf-8 -*-
# __author__ = 'shicaiping'
import os,glob
import shutil
from biocluster.core.exceptions import OptionError
from biocluster.agent import Agent
from biocluster.tool import Tool
from biocluster.config import Config
from mbio.packages.whole_transcriptome.utils import runcmd
import unittest


class EpiCovAgent(Agent):
    def __init__(self, parent):
        super(EpiCovAgent, self).__init__(parent)
        options = [
            {"name": "primer", "type": "infile", "format": "ref_rna_v2.common"},  # 引物序列文件
            {"name": "fasta", "type": "infile", "format": "ref_rna_v2.common"},  # 基因组序列文件
        ]
        self.add_option(options)
        self.split_fasta = self.add_tool("tool_lab.epicov.split_fasta")
        self.step.add_steps("epicov")
        self.on('start', self.stepstart)
        self.on('end', self.stepfinish)

    def stepstart(self):
        self.step.epicov.start()
        self.step.update()

    def stepfinish(self):
        self.step.epicov.finish()
        self.step.update()

    def check_options(self):
        """
        重写参数检测函数
        :return:
        """
        if not self.option('primer').is_set:
            raise OptionError('必须输入引物序列文件')
        if not self.option('fasta').is_set:
            raise OptionError('必须输入基因组序列文件')
        return True

    def set_resource(self):
        """
        设置所需资源，需在之类中重写此方法 self._cpu ,self._memory
        :return:
        """
        self._cpu = 1
        self._memory = "10G"

    def end(self):
        result_dir = self.add_upload_dir(self.output_dir)
        super(EpiCovAgent, self).end()


class EpiCovTool(Tool):
    def __init__(self, config):
        super(EpiCovTool, self).__init__(config)
        self.program = {
            'python': 'miniconda2/bin/python'
        }

    def run(self):
        """
        运行
        :return:
        """
        super(EpiCovTool, self).run()
        self.run_split_fasta()
        self.set_output()
        self.end()

    def run_split_fasta(self):
        options = {
            'fasta': self.option("fasta").prop['path'],
        }
        self.split_fasta.set_options(options)
        self.split_fasta.run()


    def set_output(self):
        """
        将结果文件复制到output文件夹下面
        :return:
        """
        pass


class TestFunction(unittest.TestCase):
    """
    This is test for the tool. Just run this script to do test.
    """
    def test(self):
        import random
        from mbio.workflows.single import SingleWorkflow
        from biocluster.wsheet import Sheet
        data = {
            "id": "exp_units_convert_" + str(random.randint(1, 10000)),
            "type": "tool",
            "name": "tool_lab.exp_units_convert",
            "options": dict(
                exp_matrix="/mnt/ilustre/users/sanger-dev/sg-users/shicaiping/de_tools/known_seqs_count.matrix",
                convert_type="count2tpm",
                gene_length='/mnt/ilustre/users/sanger-dev/sg-users/shicaiping/de_tools/gene_length.txt',
                float_num=4,
                intersect=True,
            )
        }
        wsheet = Sheet(data=data)
        wf = SingleWorkflow(wsheet)
        wf.run()


if __name__ == '__main__':
    unittest.main()
