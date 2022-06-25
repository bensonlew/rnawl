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


class SeqFilterAgent(Agent):
    """
    Used for fasta file filter .
    """
    def __init__(self, parent):
        super(SeqFilterAgent, self).__init__(parent)
        options = [
            {"name": "fasta_file", "type": "infile", "format": "ref_rna_v2.fasta"},  # 参考基因文件
            {"name": "max_len", "type": "int","default":100000000 },  # 过滤条件：最大基因长度
            {"name": "min_len", "type": "int","default": 0},  # 过滤条件：最小基因长度
        ]
        self.add_option(options)
        self.step.add_steps("seq_filter")
        self.on('start', self.stepstart)
        self.on('end', self.stepfinish)

    def stepstart(self):
        self.step.seq_filter.start()
        self.step.update()

    def stepfinish(self):
        self.step.seq_filter.finish()
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
        super(SeqFilterAgent, self).end()


class SeqFilterTool(Tool):
    def __init__(self, config):
        super(SeqFilterTool, self).__init__(config)
        self._version = "v1.0.1"
        self.perl =  'miniconda2/bin/'
        self.tool_path=self.config.PACKAGE_DIR+"/tool_lab/seq_filter.pl"

    def run(self):
        """
        运行
        :return:
        """
        super(SeqFilterTool, self).run()
        self.run_seq_filter()
        self.set_output()
        self.end()

    def run_seq_filter(self):
        self.logger.info("开始对fasta文件进行过滤")
        result_path = os.path.join(self.work_dir,"filtered.fasta")
        filter_cmd = '%sperl %s %s Y %d %d %s' % (self.perl, self.tool_path, self.option("fasta_file").prop["path"],self.option("min_len"),self.option("max_len"), result_path)
        command = self.add_command("seq_filter", filter_cmd, ignore_error=True)
        command.run()
        self.wait(command)
        if command.return_code == 0:
            self.logger.info("运行seq_filter完成")
        else:
            self.set_error("运行seq_filter运行出错!")
            return False



    def set_output(self):
        """
        将结果文件复制到output文件夹下面
        :return:
        """
        self.logger.info("设置结果目录")
        try:
                result=os.path.join(self.work_dir, "filtered.fasta")
                link = os.path.join(self.output_dir, "filtered.fasta")
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
            "id": "seq_filter" + str(random.randint(1, 10000)),
            "type": "tool",
            "name": "tool_lab.seq_filter",
            "instant": True,
            "options": dict(
                # fasta_file=test_dir + "/" + "dna/Saccharomyces_cerevisiae.dna.toplevel.fa",
                fasta_file ="/mnt/ilustre/users/sanger-dev/sg-users/fuwenyao/tool_lab/seq_filter/input/example.fasta" ,
                min_len=500,
                max_len=10000
            )
        }
        wsheet = Sheet(data=data)
        wf = SingleWorkflow(wsheet)
        wf.run()


if __name__ == '__main__':
    unittest.main()
