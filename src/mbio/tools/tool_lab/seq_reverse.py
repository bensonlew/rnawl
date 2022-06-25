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


class SeqReverseAgent(Agent):
    """
    Used for fasta file seq reverse .

    detail:
        remethod-r:raw:ATTC new:CTTA
        remethod-c:raw:ATTC new:TAAG
        remethod-cr:raw:ATTC new:GAAT

    """
    def __init__(self, parent):
        super(SeqReverseAgent, self).__init__(parent)
        options = [
            {"name": "fasta_file", "type": "infile", "format": "ref_rna_v2.fasta"},  # 参考基因文件
            {"name": "re_method", "type": "string","default":"rc" },  # 序列转化的方法，默认c，即直接取["c","r","cr"]
        ]
        self.add_option(options)
        self.step.add_steps("seq_reverse")
        self.on('start', self.stepstart)
        self.on('end', self.stepfinish)

    def stepstart(self):
        self.step.seq_reverse.start()
        self.step.update()

    def stepfinish(self):
        self.step.seq_reverse.finish()
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
        super(SeqReverseAgent, self).end()


class SeqReverseTool(Tool):
    def __init__(self, config):
        super(SeqReverseTool, self).__init__(config)
        self._version = "v1.0.1"
        self.python_path = '/miniconda2/bin/python'
        self.perl =  '/miniconda2/bin/perl'
        self.tool_path_old=self.config.PACKAGE_DIR+"/tool_lab/seq_reverse.py"
        self.tool_path= self.config.PACKAGE_DIR + "/tool_lab/seq_reverse.pl"

    def run(self):
        """
        运行
        :return:
        """
        super(SeqReverseTool, self).run()
        self.run_seq_reverse()
        self.set_output()
        self.end()

    def run_seq_reverse(self):
        self.logger.info("开始对fasta文件进行序列取反")
        reversed_cmd = '%s %s -input %s -type %s -output %s ' % (self.perl, self.tool_path, self.option("fasta_file").prop["path"],self.option("re_method"),self.output_dir)
        command = self.add_command("seq_reverse", reversed_cmd, ignore_error=True)
        command.run()
        self.wait(command)
        if command.return_code == 0:
            self.logger.info("运行reversed_cmd完成")
        else:
            self.set_error("运行reversed_cmd运行出错!")
            return False



    def set_output(self):
        """
        将结果文件复制到output文件夹下面
        :return:
        """
        # self.logger.info("设置结果目录")
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
            "id": "seq_reverse" + str(random.randint(1, 10000)),
            "type": "tool",
            "name": "tool_lab.seq_reverse",
            "instant": True,
            "options": dict(
                # fasta_file=test_dir + "/" + "dna/Saccharomyces_cerevisiae.dna.toplevel.fa",
                fasta_file ="/mnt/ilustre/users/sanger-dev/sg-users/fuwenyao/tool_lab/seq_reverse/test2.fasta" ,
                # min_len=500,
                # max_len=10000
            )
        }
        wsheet = Sheet(data=data)
        wf = SingleWorkflow(wsheet)
        wf.run()


if __name__ == '__main__':
    unittest.main()
