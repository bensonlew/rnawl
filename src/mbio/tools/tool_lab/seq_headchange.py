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


class SeqHeadchangeAgent(Agent):
    """
    Used for  change fasta header
    input:
    >seq1
    ATCGTCGATCGA
    >seq1
    TCGAATCGTCATCGTC

    symbol:"seq1"  connect:"_"   position:"back"
    out:
    >S2_seq1
    ATCGTCGATCGA
    >S2_seq1
    TCGAATCGTCATCGTC

    """
    def __init__(self, parent):
        super(SeqHeadchangeAgent, self).__init__(parent)
        options = [
            # 参考基因文件
            {"name": "fasta_file", "type": "infile", "format": "ref_rna_v2.fasta"},
            #标识符添加位置  ["back","head"]
            {"name": "position", "type": "string","default":"before" },
            #连接符
            {"name": "connect", "type": "string", "default": "_"},
            # 标识符
            {"name": "symbol", "type": "string"},

        ]
        self.add_option(options)
        self.step.add_steps("seq_headchange")
        self.on('start', self.stepstart)
        self.on('end', self.stepfinish)

    def stepstart(self):
        self.step.seq_headchange.start()
        self.step.update()

    def stepfinish(self):
        self.step.seq_headchange.finish()
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
        super(SeqHeadchangeAgent, self).end()


class SeqHeadchangeTool(Tool):
    def __init__(self, config):
        super(SeqHeadchangeTool, self).__init__(config)
        self._version = "v1.0.1"
        self.python_path = '/program/Python/bin/python'
        # self.perl =  'program/perl/perls/perl-5.24.0/bin/'
        self.perl = self.config.SOFTWARE_DIR+'/program/perl/perls/perl-5.24.0/bin/'
        self.tool_path=self.config.PACKAGE_DIR+"/tool_lab/seq_headchange.pl"

    def run(self):
        """
        运行
        :return:
        """
        super(SeqHeadchangeTool, self).run()
        self.run_seq_headchange()
        self.set_output()
        self.end()

    def run_seq_headchange(self):
        self.logger.info("开始对fasta文件进行序列取反")
        outfa_dir=os.path.join(self.work_dir,"header_changed.fasta")
        run_seq_headchange_cmd = '{}perl {} -fa {} -symbol {} -position {} -connect "{}" -out {}'\
                       .format(self.perl, self.tool_path, self.option("fasta_file").prop["path"],
                               self.option("symbol"),self.option("position"),self.option("connect"),outfa_dir)
        command = self.add_command("seq_headchange", run_seq_headchange_cmd, ignore_error=True,shell=True)
        command.run()
        self.wait(command)
        if command.return_code == 0:
            self.logger.info("运行seq_headchange完成")
        else:
            self.set_error("运行seq_headchange运行出错!")
            return False



    def set_output(self):
        """
        将结果文件复制到output文件夹下面
        :return:
        """
        self.logger.info("设置结果目录")
        try:
                result=os.path.join(self.work_dir, "header_changed.fasta")
                link = os.path.join(self.output_dir, "header_changed.fasta")
                if os.path.exists(link):
                    os.remove(link)
                os.link(result, link)
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
        test_dir = "/mnt/ilustre/users/sanger-dev/app/database/Genome_DB_finish/fungi/Saccharomyces_cerevisiae/Ensemble_release_39/"
        data = {
            "id": "seq_headchange" + str(random.randint(1, 10000)),
            "type": "tool",
            "name": "tool_lab.seq_headchange",
            "instant": True,
            "options": dict(
                # fasta_file=test_dir + "/" + "dna/Saccharomyces_cerevisiae.dna.toplevel.fa",
                fasta_file ="/mnt/ilustre/users/sanger-dev/sg-users/fuwenyao/tool_lab/FAheader/input/ref.fa" ,
                # table="Standard",
                position="before",
                symbol = "pppp",
                connect="_"

                # min_len=500,
                # max_len=10000
            )
        }
        wsheet = Sheet(data=data)
        wf = SingleWorkflow(wsheet)
        wf.run()


if __name__ == '__main__':
    unittest.main()
