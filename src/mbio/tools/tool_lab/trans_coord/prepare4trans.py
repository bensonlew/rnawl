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


class Prepare4transAgent(Agent):
    """
    该tool

    """
    def __init__(self, parent):
        super(Prepare4transAgent, self).__init__(parent)
        options = [
            {"name": "raw_fasta_file", "type": "infile", "format": "ref_rna_v2.common"},  # 参考基因文件
            {"name": "target_fasta_file", "type": "infile", "format": "ref_rna_v2.common"},  # 参考基因文件
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
        if not self.option('raw_fasta_file').is_set:
            raise OptionError('必须输入当前参考基因组序列文件')
        if not self.option('target_fasta_file').is_set:
            raise OptionError('必须输入当前参考基因组序列文件')
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
        super(Prepare4transAgent, self).end()


class Prepare4transTool(Tool):
    def __init__(self, config):
        super(Prepare4transTool, self).__init__(config)
        self._version = "v1.0.1"
        self.python_path = '/miniconda2/bin/python'
        self.perl =  'miniconda2/bin/'
        self.tool_path=self.config.PACKAGE_DIR+"/tool_lab/trans_coord/split_fasta_sgene.py"
        self.ucsctool_path =  "/bioinfo/align/ucsc_tools/"
        self.lastz_path = self.config.SOFTWARE_DIR + "/miniconda2/bin"
        # self.parafly = "/program/parafly-r2013-01-21/src/ParaFly"
        self.set_environ(PATH=self.ucsctool_path)
        self.set_environ(PATH=self.lastz_path)

    def run(self):
        """
        运行
        :return:
        """
        super(Prepare4transTool, self).run()
        self.run_target_fasplit()
        self.run_query_fasplit()
        self.run_tfaToTwoBit()
        self.run_qfaToTwoBit()
        self.run_ttwoBitInfo()
        self.run_qtwoBitInfo()
        self.set_output()
        self.end()

    def run_target_fasplit(self):
        # Info_dir = os.path.join(self.output_dir, "Info")
        # # lastz = os.path.join(self.output_dir, "lastz")
        # os.mkdir(Info_dir)
        tfasplit_cmd = '{}faSplit -lift={}/target.lft size {} ' \
                      '-oneFile 5000000  -extra=10000 {}/target'.format(self.ucsctool_path, self.output_dir,
                                                                             self.option("raw_fasta_file").prop["path"], self.output_dir)
        command = self.add_command("tfasplit", tfasplit_cmd, ignore_error=True)
        command.run()
        self.wait(command)
        if command.return_code == 0:
            self.logger.info("运行fasplitsplit_cmd完成")
        else:
            self.set_error("运行fasplit运行出错!")
            return False

    def run_query_fasplit(self):
        qfasplit_cmd = '{}faSplit -lift={}/query.lft size {} ' \
                      '-oneFile 5000000  -extra=10000 {}/query'.format(self.ucsctool_path, self.output_dir,
                                                                             self.option("target_fasta_file").prop["path"],
                                                                             self.output_dir)
        command = self.add_command("qfasplit", qfasplit_cmd, ignore_error=True)
        command.run()
        self.wait(command)
        if command.return_code == 0:
            self.logger.info("运行qfasplit完成")
        else:
            self.set_error("qfasplit运行出错!")
            return False

    def run_tfaToTwoBit(self):
        tfatotwobit_cmd='{}faToTwoBit {}  {}/target.2bit'.format(self.ucsctool_path,self.option("raw_fasta_file").prop["path"],self.output_dir)
        command = self.add_command("tfatotwobit_cmd", tfatotwobit_cmd, ignore_error=True)
        command.run()
        self.wait(command)
        if command.return_code == 0:
            self.logger.info("运行tfatotwobit_cmd完成")
        else:
            self.set_error("tfatotwobit_cmd运行出错!")
            return False

    def run_qfaToTwoBit(self):
        qfatotwobit_cmd = '{}faToTwoBit {}  {}/query.2bit'.format(self.ucsctool_path,
                                                                        self.option("target_fasta_file").prop["path"],
                                                                        self.output_dir)
        command = self.add_command("qfatotwobit_cmd", qfatotwobit_cmd, ignore_error=True)
        command.run()
        self.wait(command)
        if command.return_code == 0:
            self.logger.info("运行qfatotwobit_cmd完成")
        else:
            self.set_error("qfatotwobit_cmd运行出错!")
            return False

    def run_ttwoBitInfo(self):
        ttwoBitInfo_cmd='{}twoBitInfo {}/target.2bit {}/target.chrom.sizes'.format(self.ucsctool_path,self.output_dir,self.output_dir)
        command = self.add_command("ttwobitinfo_cmd", ttwoBitInfo_cmd, ignore_error=True)
        command.run()
        self.wait(command)
        if command.return_code == 0:
            self.logger.info("运行ttwoBitInfo_cmd完成")
        else:
            self.set_error("ttwoBitInfo_cmd运行出错!")
            return False

    def run_qtwoBitInfo(self):
        qtwoBitInfo_cmd = '{}twoBitInfo {}/query.2bit {}/query.chrom.sizes'.format(self.ucsctool_path,
                                                                                              self.output_dir,
                                                                                              self.output_dir)
        command = self.add_command("qtwobitinfo_cmd", qtwoBitInfo_cmd, ignore_error=True)
        command.run()
        self.wait(command)
        if command.return_code == 0:
            self.logger.info("运行qtwoBitInfo_cmd完成")
        else:
            self.set_error("qtwoBitInfo_cmd运行出错!")
            return False



    def set_output(self):
        """
        将结果文件复制到output文件夹下面
        :return:
        """
        self.logger.info("设置结果目录")



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
            "name": "tool_lab.trans_coord.prepare4trans",
            "instant": True,
            "options": dict(
                # fasta_file=test_dir + "/" + "dna/Saccharomyces_cerevisiae.dna.toplevel.fa",
                raw_fasta_file ="/mnt/ilustre/users/sanger-dev/sg-users/fuwenyao/tool_lab/TransCoord/input/ASM14920v2.fasta" ,
                target_fasta_file="/mnt/ilustre/users/sanger-dev/sg-users/fuwenyao/tool_lab/TransCoord/input/ASM18445v3.fasta"
                # min_len=500,
                # max_len=10000
            )
        }
        wsheet = Sheet(data=data)
        wf = SingleWorkflow(wsheet)
        wf.run()


if __name__ == '__main__':
    unittest.main()
