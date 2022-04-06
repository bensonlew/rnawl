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


class OrffinderAgent(Agent):
    """
    在目标DNA序列中搜寻开发阅读框，可以输出每个ORF所在的区域，并翻译成对应的蛋白序列
    此工具可以为新预测的DNA序列查找潜在的蛋白编码区
    """
    def __init__(self, parent):
        super(OrffinderAgent, self).__init__(parent)
        options = [
            {"name": "input_file", "type": "infile", "format": "sequence.fasta"},  # FASTA序列文件
            {"name": "start_codon", "type": "int", "default": 2},
            # 0 = "ATG" only
            # 1 = "ATG" and alternative initiation codons
            # 2 = any sense codon
            {"name": "start_pos", "type": "string", "default": '1,2,3'},
            {"name": "genetic_code", "type": "int", "default": 1},
            # Genetic code to use (1-31)
            # see https://www.ncbi.nlm.nih.gov/Taxonomy/Utils/wprintgc.cgi for details
            {"name": "minimal_length", "type": "int", "default": 75},
            # Minimal length of the ORF (nt)
            # Value less than 30 is automatically changed by 30.
            {"name": "strand", "type": "string", "default": "both"},
            # Output ORFs on specified strand only (both|plus|minus)
        ]
        self.add_option(options)
        self.step.add_steps("orffinder")
        self.on('start', self.stepstart)
        self.on('end', self.stepfinish)

    def stepstart(self):
        self.step.orffinder.start()
        self.step.update()

    def stepfinish(self):
        self.step.orffinder.finish()
        self.step.update()

    def check_options(self):
        if not self.option("input_file").is_set:
            raise OptionError('请输入FASTA文件')
        if self.option("start_codon") not in [0, 1, 2]:
            raise OptionError('开放阅读框参数输入错误')
        if self.option("strand") not in ['minus', 'plus', 'both']:
            raise OptionError('链参数输入错误')
        for start_pos in self.option("start_pos").split(","):
            if start_pos not in ['1', '2', '3']:
                raise OptionError('查找ORF阅读框参数输入错误')
        if self.option("genetic_code") not in range(1, 32):
            raise OptionError('密码子表参数输入错误')
        return True

    def set_resource(self):
        """
        设置所需资源，需在之类中重写此方法 self._cpu ,self._memory
        :return:
        """
        self._cpu = 1
        self._memory = "5G"

    def end(self):
        result_dir = self.add_upload_dir(self.output_dir)
        super(OrffinderAgent, self).end()


class OrffinderTool(Tool):
    def __init__(self, config):
        super(OrffinderTool, self).__init__(config)
        self.program = {
            'python': 'miniconda2/bin/python',
            'orffinder': 'bioinfo/rna/ORFfinder'
        }


    def run(self):
        """
        运行
        :return:
        """
        super(OrffinderTool, self).run()
        for start_pos in self.option("start_pos").split(","):
            self.run_orffinder(int(start_pos))
        # self.set_output()
        self.end()

    def run_orffinder(self, start_pos):
        out_file = os.path.join(self.output_dir, "orf_" + self.option('strand') + "_" + str(start_pos) + ".fa")
        cmd = '{}'.format(self.program['orffinder'])
        cmd += ' -in {}'.format(self.option('input_file').prop['path'])
        cmd += ' -s {}'.format(self.option('start_codon'))
        cmd += ' -ml {}'.format(self.option('minimal_length'))
        cmd += ' -strand {}'.format(self.option('strand'))
        cmd += ' -b {}'.format(start_pos)
        cmd += ' -outfmt 1 -out {}'.format(out_file)
        cmd_name = 'orffinder_' + str(start_pos)
        runcmd(self, cmd_name, cmd)

    def set_output(self):
        """
        将结果文件复制到output文件夹下面
        :return:
        """
        self.logger.info("设置结果目录")
        try:
            output_files = glob.glob("page*")
            for file in output_files:
                os.link(os.path.join(self.work_dir, file), os.path.join(self.output_dir, file))
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
        data = {
            "id": "orffinder_" + str(random.randint(1, 10000)),
            "type": "tool",
            "name": "tool_lab.orffinder",
            "options": dict(
                input_file="/mnt/lustre/users/sanger/sg-users/shicaiping/FASTA_example.fsa",
                start_codon=2,
                start_pos='1,2,3',
                genetic_code=1,
                minimal_length=75,
                strand="both",
            )
        }
        wsheet = Sheet(data=data)
        wf = SingleWorkflow(wsheet)
        wf.run()


if __name__ == '__main__':
    unittest.main()
