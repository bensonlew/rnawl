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


class Fa2chainAgent(Agent):
    """
    该tool用于生成切割后的fasta文件与queryfasta相关的chain文件

    """
    def __init__(self, parent):
        super(Fa2chainAgent, self).__init__(parent)
        options = [
            {"name": "fasta_file", "type": "infile", "format": "ref_rna_v2.common"},  #
            {"name": "info_dir", "type": "infile", "format": "ref_rna_v2.common_dir"},  #
        ]
        self.add_option(options)
        self.step.add_steps("fa2chain")
        self.on('start', self.stepstart)
        self.on('end', self.stepfinish)

    def stepstart(self):
        self.step.fa2chain.start()
        self.step.update()

    def stepfinish(self):
        self.step.fa2chain.finish()
        self.step.update()

    def check_options(self):
        """
        重写参数检测函数
        :return:
        """
        if not self.option('fasta_file').is_set:
            raise OptionError('必须输入分割后的fasta文件')
        if not self.option('info_dir').is_set:
            raise OptionError('必须输入当前info文件夹')
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
        super(Fa2chainAgent, self).end()


class Fa2chainTool(Tool):
    def __init__(self, config):
        super(Fa2chainTool, self).__init__(config)
        self._version = "v1.0.1"
        self.python_path = '/miniconda2/bin/python'
        self.perl =  'program/perl/perls/perl-5.24.0/bin/'
        self.tool_path=self.config.PACKAGE_DIR+"/tool_lab/trans_coord/split_fasta_sgene.py"
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
        super(Fa2chainTool, self).run()
        self.run_lastz()
        self.run_lavToPsl()
        self.run_liftup()
        self.run_axtChain()
        self.set_output()
        self.end()

    def run_lastz(self):
        self.lastz_dir = os.path.join(self.output_dir, "lastz")
        os.mkdir(self.lastz_dir)
        self.ch_name=os.path.basename(self.option("fasta_file").prop["path"])
        lastz_cmd = '{}lastz {} {}/query.fa --masking=50' \
                    ' --hspthresh=2200 --ydrop=3400 --gappedthresh=4000 --inner=2000 --ambiguous=iupac ' \
                    '--output={}/{}.query.lav'\
            .format(self.lastz_path, self.option('fasta_file').prop["path"],self.option("info_dir").prop["path"],self.lastz_dir,self.ch_name)
        command = self.add_command("lastz_cmd", lastz_cmd, ignore_error=True)
        command.run()
        self.wait(command)
        if command.return_code == 0:
            self.logger.info("运行lastz_cmd完成")
        else:
            self.set_error("运行lastz_cmd运行出错!")
            return False

    def run_lavToPsl(self):
        self.psl_dir=os.path.join(self.output_dir,"psl")
        self.psl_detail_dir=os.path.join(self.output_dir,"psl/psl")
        self.chain_dir=os.path.join(self.output_dir,"chain")
        # self.chain_detail_dir=os.path.join(self.output_dir,"chain/chain")
        os.mkdir(self.psl_dir)
        os.mkdir(self.psl_detail_dir)
        os.mkdir(self.chain_dir)
        # os.mkdir(self.chain_detail_dir)
        lavtopsl_cmd = "{}lavToPsl {}/{}.query.lav  {}/{}.query.lav.psl"\
            .format(self.ucsctool_path,self.lastz_dir,self.ch_name,self.psl_dir,self.ch_name)
        command = self.add_command("lavtopsl_cmd", lavtopsl_cmd, ignore_error=True)
        command.run()
        self.wait(command)
        if command.return_code == 0:
            self.logger.info("运行lavtopsl_cmd完成")
        else:
            self.set_error("运行lavtopsl_cmd运行出错!")
            return False

    def run_liftup(self):
        liftup_cmd="{}liftUp -pslQ {}/query.{}.query.lav.psl {}/query.lft " \
                   "warn {}/{}.query.lav.psl".format(self.ucsctool_path ,
                                                     self.psl_detail_dir,self.ch_name,self.option("info_dir").prop["path"],self.psl_dir,self.ch_name )
        command = self.add_command("liftup_cmd", liftup_cmd, ignore_error=True)
        command.run()
        self.wait(command)
        if command.return_code == 0:
            self.logger.info("运行liftup_cmd完成")
        else:
            self.set_error("运行liftup_cmd运行出错!")
            return False

    def run_axtChain(self):
        axtchain_cmd="{}axtChain -linearGap=medium -psl {}/query.{}.query.lav.psl " \
                     "{}/target.2bit {}/query.2bit " \
                     "{}/query.{}.query.lav.psl.chain"\
            .format(self.ucsctool_path,self.psl_detail_dir,self.ch_name,self.option("info_dir").prop["path"],
                    self.option("info_dir").prop["path"],self.chain_dir,self.ch_name)
        command = self.add_command("axtchain_cmd", axtchain_cmd, ignore_error=True)
        command.run()
        self.wait(command)
        if command.return_code == 0:
            self.logger.info("运行axtchain_cmd完成")
        else:
            self.set_error("运行axtchain_cmd运行出错!")
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
            "id": "fa2chain" + str(random.randint(1, 10000)),
            "type": "tool",
            "name": "tool_lab.trans_coord.fa2chain",
            "instant": True,
            "options": dict(
                # fasta_file=test_dir + "/" + "dna/Saccharomyces_cerevisiae.dna.toplevel.fa",
                fasta_file ="/mnt/ilustre/users/sanger-dev/workspace/20200426/Single_trans_coord3482/TransCoord/fasta_split/chr_9.fa" ,
                info_dir="/mnt/ilustre/users/sanger-dev/workspace/20200426/Single_trans_coord3482/TransCoord/Prepare4trans/output/"
                # min_len=500,
                # max_len=10000
            )
        }
        wsheet = Sheet(data=data)
        wf = SingleWorkflow(wsheet)
        wf.run()


if __name__ == '__main__':
    unittest.main()
