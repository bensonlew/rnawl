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


class NetchainsubsetAgent(Agent):
    """
    net文件的合并

    """
    def __init__(self, parent):
        super(NetchainsubsetAgent, self).__init__(parent)
        options = [
            {"name": "net_dir", "type": "infile", "format": "ref_rna_v2.common_dir"},  #
            {"name": "target_chain", "type": "infile", "format": "ref_rna_v2.common"},  #
            # {"name": "info_dir", "type": "infile", "format": "ref_rna_v2.common_dir"},  #
        ]
        self.add_option(options)
        self.step.add_steps("net_chain")
        self.on('start', self.stepstart)
        self.on('end', self.stepfinish)

    def stepstart(self):
        self.step.net_chain.start()
        self.step.update()

    def stepfinish(self):
        self.step.net_chain.finish()
        self.step.update()

    def check_options(self):
        """
        重写参数检测函数
        :return:
        """
        if not self.option('net_dir').is_set:
            raise OptionError('必须输入net_dir')
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
        super(NetchainsubsetAgent, self).end()


class NetchainsubsetTool(Tool):
    def __init__(self, config):
        super(NetchainsubsetTool, self).__init__(config)
        self._version = "v1.0.1"
        self.python_path = '/miniconda2/bin/python'
        self.perl =  'program/perl/perls/perl-5.24.0/bin/'
        self.tool_path=self.config.PACKAGE_DIR+"/tool_lab/trans_coord/split_fasta_sgene.py"
        self.chainMergeSort_path = self.config.SOFTWARE_DIR+ "/bioinfo/align/ucsc_tools/"
        self.ucsctool_path =  "/bioinfo/align/ucsc_tools/"
        self.lastz_path = "/bioinfo/ref_rna_v2/miniconda2/bin/"
        # self.parafly = "/program/parafly-r2013-01-21/src/ParaFly"
        self.set_environ(PATH=self.chainMergeSort_path)
        self.set_environ(PATH=self.lastz_path)

    def run(self):
        """
        运行
        :return:
        """
        super(NetchainsubsetTool, self).run()
        self.cat_nets()
        self.format_adjust()
        self.set_output()
        self.end()

    def cat_nets(self):
        subprocess.call("cat {}/*.net > target.net".format(self.option('net_dir').prop["path"]),shell=True)
        self.run_netChainSubset()

    def run_netChainSubset(self):
        netchain_cmd= '{}netChainSubset target.net {} target.chain'\
            .format(self.ucsctool_path,self.option("target_chain").prop["path"])
        command = self.add_command("netchain_cmd", netchain_cmd, ignore_error=True)
        command.run()
        self.wait(command)
        if command.return_code == 0:
            self.logger.info("运行netchain_cmd完成")
        else:
            self.set_error("运行netchain_cmd运行出错!")
            return False

    def format_adjust(self):
        jobs="less target.chain |grep -v \"#\" > {}/lift.chain".format(self.output_dir)
        subprocess.call(jobs,shell=True)

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
            "id": "netchainsubset" + str(random.randint(1, 10000)),
            "type": "tool",
            "name": "tool_lab.trans_coord.netchainsubset",
            "instant": True,
            "options": dict(
                # fasta_file=test_dir + "/" + "dna/Saccharomyces_cerevisiae.dna.toplevel.fa",
                net_dir ="/mnt/ilustre/users/sanger-dev/sg-users/fuwenyao/tool_lab/TransCoord/output1/Step1.chain/net" ,
                target_chain="/mnt/ilustre/users/sanger-dev/sg-users/fuwenyao/tool_lab/TransCoord/output1/Step1.chain/chain/chain/target.chain"
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
