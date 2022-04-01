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


class SeqStatAgent(Agent):
    """
    Used for fasta seq stat .
    """
    def __init__(self, parent):
        super(SeqStatAgent, self).__init__(parent)
        options = [
            {"name": "fasta_file", "type": "infile", "format": "ref_rna_v2.fasta"},  # 参考基因组文件
            # {"name": "step", "type": "int","default": 500 },  # 统计步长    统计的步长
            {"name": "group_num", "type": "int", "default":10},  # 统计组数
            {"name": "min_len", "type": "int", "default":0},  # 最小统计长度   该长度一下的reads不进行统计
        ]
        self.add_option(options)
        self.step.add_steps("seq_stat")
        self.on('start', self.stepstart)
        self.on('end', self.stepfinish)

    def stepstart(self):
        self.step.seq_stat.start()
        self.step.update()

    def stepfinish(self):
        self.step.seq_stat.finish()
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
        super(SeqStatAgent, self).end()


class SeqStatTool(Tool):
    def __init__(self, config):
        super(SeqStatTool, self).__init__(config)
        self._version = "v1.0.1"
        self.python_path = '/program/Python/bin/python'
        self.perl = 'program/perl/perls/perl-5.24.0/bin/perl'
        self.tool_path=self.config.PACKAGE_DIR+"/tool_lab/seq_stat.pl"

    def run(self):
        """
        运行
        :return:
        """
        super(SeqStatTool, self).run()
        self.run_seq_stat()
        self.set_output()
        self.end()

    def run_seq_stat(self):
        self.logger.info("开始对fasta文件进行统计")
        stat_cmd = '{} {} '.format(self.perl,self.tool_path)
        stat_cmd += '-input {} '.format(self.option("fasta_file").prop["path"])
        stat_cmd += '-split {} '.format(self.option("group_num"))
        # stat_cmd += '-s {} '.format(self.option("step"))
        stat_cmd += '-output seq_stat '
        command = self.add_command("seq_stat", stat_cmd, ignore_error=True)
        command.run()
        self.wait(command)
        if command.return_code == 0:
            self.logger.info("运行seq_stat完成")
        else:
            self.set_error("运行seq_stat运行出错!")
            return False



    def set_output(self):
        """
        将结果文件复制到output文件夹下面
        :return:
        """
        self.logger.info("设置结果目录")
        try:
                result_detail=os.path.join(self.work_dir, "seq_stat.length")
                result_stat = os.path.join(self.work_dir, "seq_stat.distribution")
                link_detail = os.path.join(self.output_dir, "seq_stat.length")
                link_stat = os.path.join(self.output_dir, "seq_stat.distribution")
                if os.path.exists(link_detail):
                    os.remove(link_detail)
                os.link(result_detail, link_detail)
                if os.path.exists(link_stat):
                    os.remove(link_stat)
                os.link(result_stat, link_stat)
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
            "id": "seq_stat" + str(random.randint(1, 10000)),
            "type": "tool",
            "name": "tool_lab.seq_stat",
            "instant": True,
            "options": dict(
                # fasta_file=test_dir + "/" + "dna/Saccharomyces_cerevisiae.dna.toplevel.fa",
                fasta_file ="/mnt/ilustre/users/sanger-dev/sg-users/fuwenyao/tool_lab/seq_statistic_new/input/example.fasta" ,
                group_num=10,
                # step=500
            )
        }
        wsheet = Sheet(data=data)
        wf = SingleWorkflow(wsheet)
        wf.run()


if __name__ == '__main__':
    unittest.main()
