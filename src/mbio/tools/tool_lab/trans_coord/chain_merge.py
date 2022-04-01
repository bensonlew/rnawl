# -*- coding: utf-8 -*-
# __author__ = 'fwy'
import os
import shutil
from biocluster.core.exceptions import OptionError
from biocluster.agent import Agent
from biocluster.tool import Tool
from biocluster.config import Config
import subprocess
import re
import pandas as pd
import unittest


class ChainMergeAgent(Agent):
    """
    多余所有得到的chain文件进行排序分割和合并

    """
    def __init__(self, parent):
        super(ChainMergeAgent, self).__init__(parent)
        options = [
            {"name": "chain_dir", "type": "infile", "format": "ref_rna_v2.common_dir"},  #
            # {"name": "info_dir", "type": "infile", "format": "ref_rna_v2.common_dir"},  #
        ]
        self.add_option(options)
        self.step.add_steps("chain_merge")
        self.on('start', self.stepstart)
        self.on('end', self.stepfinish)

    def stepstart(self):
        self.step.chain_merge.start()
        self.step.update()

    def stepfinish(self):
        self.step.chain_merge.finish()
        self.step.update()

    def check_options(self):
        """
        重写参数检测函数
        :return:
        """
        if not self.option('chain_dir').is_set:
            raise OptionError('必须输入chain_dir')
        return True

    def set_resource(self):
        """
        设置所需资源，需在之类中重写此方法 self._cpu ,self._memory
        :return:
        """
        self._cpu = 2
        self._memory = "10G"

    def end(self):
        result_dir = self.add_upload_dir(self.output_dir)
        super(ChainMergeAgent, self).end()


class ChainMergeTool(Tool):
    def __init__(self, config):
        super(ChainMergeTool, self).__init__(config)
        self._version = "v1.0.1"
        self.python_path = '/program/Python/bin/python'
        self.perl =  'program/perl/perls/perl-5.24.0/bin/'
        self.tool_path=self.config.PACKAGE_DIR+"/tool_lab/trans_coord/split_fasta_sgene.py"
        self.chainMergeSort_path = self.config.SOFTWARE_DIR+ "/bioinfo/align/ucsc_tools/"
        self.ucsctool_path =  "/bioinfo/align/ucsc_tools/"
        self.lastz_path = "/bioinfo/ref_rna_v2/miniconda2/bin/"
        # self.parafly = "/program/parafly-r2013-01-21/src/ParaFly"
        self.set_environ(PATH=self.ucsctool_path)
        self.set_environ(PATH=self.lastz_path)

    def run(self):
        """
        运行
        :return:
        """
        super(ChainMergeTool, self).run()
        self.run_chain_merge()
        self.run_chainSplit()
        self.final_chain_marge()
        self.set_output()
        self.end()

    def run_chain_merge(self):
        os.popen("{}chainMergeSort {}/*.chain > mergesort.chain".format(self.chainMergeSort_path,self.option('chain_dir').prop["path"]))
    # def run_chain_merge(self):
    #     chain_files=[]
    #     chain_list=os.listdir(self.option('chain_dir').prop["path"])
    #     for  i in chain_list:
    #         chain_files.append(os.path.join(self.option('chain_dir').prop["path"],i))
    #     pack_chain=" ".join(chain_files)
    #     chainmergesort_cmd = '{}chainMergeSort {} '.format(self.ucsctool_path,pack_chain)
    #     command = self.add_command("chainmergesort_cmd", chainmergesort_cmd, ignore_error=True)
    #     command.run()
    #     self.wait(command)
    #     if command.return_code == 0:
    #         self.logger.info("运行chainmergesort_cmd完成")
    #     else:
    #         self.set_error("运行chainmergesort_cmd运行出错!")
    #         return False

    def run_chainSplit(self):
        self.chain_dir = os.path.join(self.output_dir, "chain")
        chainsplit_cmd = '{}chainSplit {} mergesort.chain'.format(self.ucsctool_path,self.chain_dir)
        command = self.add_command("chainsplit_cmd", chainsplit_cmd, ignore_error=True)
        command.run()
        self.wait(command)
        if command.return_code == 0:
            self.logger.info("运行chainsplit_cmd完成")
        else:
            self.set_error("运行chainsplit_cmd运行出错!")
            return False

    def final_chain_marge(self):
        os.system("cat {}/*.chain > {}/target.chain".format(self.chain_dir,self.output_dir))


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
            "id": "chain_merge" + str(random.randint(1, 10000)),
            "type": "tool",
            "name": "tool_lab.trans_coord.chain_merge",
            "instant": True,
            "options": dict(
                # fasta_file=test_dir + "/" + "dna/Saccharomyces_cerevisiae.dna.toplevel.fa",
                chain_dir ="/mnt/ilustre/users/sanger-dev/sg-users/fuwenyao/tool_lab/TransCoord/output1/Step1.chain/chain" ,
                # info_dir="/mnt/ilustre/users/sanger-dev/sg-users/fuwenyao/tool_lab/TransCoord/output1/Step1.chain/Info"
                # min_len=500,
                # max_len=10000
            )
        }
        wsheet = Sheet(data=data)
        wf = SingleWorkflow(wsheet)
        wf.run()


if __name__ == '__main__':
    unittest.main()
