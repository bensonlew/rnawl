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


class FastaSplitSgeneAgent(Agent):
    """
    该tool

    """
    def __init__(self, parent):
        super(FastaSplitSgeneAgent, self).__init__(parent)
        options = [
            {"name": "fasta_file", "type": "infile", "format": "ref_rna_v2.fasta"},  # 参考基因文件
        ]
        self.add_option(options)
        self.step.add_steps("fasta_split")
        self.on('start', self.stepstart)
        self.on('end', self.stepfinish)

    def stepstart(self):
        self.step.fasta_split.start()
        self.step.update()

    def stepfinish(self):
        self.step.fasta_split.finish()
        self.step.update()

    def check_options(self):
        """
        重写参数检测函数
        :return:
        """
        if not self.option('fasta_file').is_set:
            raise OptionError('必须输入参考基因组序列文件')
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
        super(FastaSplitSgeneAgent, self).end()


class FastaSplitSgeneTool(Tool):
    def __init__(self, config):
        super(FastaSplitSgeneTool, self).__init__(config)
        self._version = "v1.0.1"
        self.python_path = '/miniconda2/bin/python'
        self.perl =  'program/perl/perls/perl-5.24.0/bin/'
        self.tool_path=self.config.PACKAGE_DIR+"/tool_lab/trans_coord/split_fasta_sgene.py"

    def run(self):
        """
        运行
        :return:
        """
        super(FastaSplitSgeneTool, self).run()
        self.run_fasta_split()
        self.set_output()
        self.end()

    def run_fasta_split(self):
        self.logger.info("开始对fasta文件进行分基因切割")
        split_cmd = '%s %s -i %s -o %s' % (self.python_path, self.tool_path, self.option("fasta_file").prop["path"],self.output_dir)
        command = self.add_command("fasta_split", split_cmd, ignore_error=True)
        command.run()
        self.wait(command)
        if command.return_code == 0:
            self.logger.info("运行split_cmd完成")
        else:
            self.set_error("运行split_cmd运行出错!")
            return False



    def set_output(self):
        """
        将结果文件复制到output文件夹下面
        :return:
        """
        self.logger.info("设置结果目录")
        # try:
        #         result=os.path.join(self.work_dir, "reversed.{}.fasta".format(self.option("re_method")))
        #         link = os.path.join(self.output_dir, "reversed.fasta")
        #         if os.path.exists(link):
        #             os.remove(link)
        #         os.link(result, link)
        # except Exception as e:
        #     self.logger.info("设置结果目录失败{}".format(e))


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
            "id": "fasta_split_sgene" + str(random.randint(1, 10000)),
            "type": "tool",
            "name": "tool_lab.trans_coord.fasta_split_sgene",
            "instant": True,
            "options": dict(
                # fasta_file=test_dir + "/" + "dna/Saccharomyces_cerevisiae.dna.toplevel.fa",
                fasta_file ="/mnt/ilustre/users/sanger-dev/sg-users/fuwenyao/tool_lab/TransCoord/input/ASM14920v2.fasta" ,
                # min_len=500,
                # max_len=10000
            )
        }
        wsheet = Sheet(data=data)
        wf = SingleWorkflow(wsheet)
        wf.run()


if __name__ == '__main__':
    unittest.main()
