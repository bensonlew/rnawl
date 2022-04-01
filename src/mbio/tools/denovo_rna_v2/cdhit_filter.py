# -*- coding: utf-8 -*-
# __author__ = 'liubinxu'
from biocluster.agent import Agent
from biocluster.tool import Tool
import os
import subprocess
from biocluster.core.exceptions import OptionError
import shutil
import unittest

class CdhitFilterAgent(Agent):
    """
    过滤转录本冗余序列
    """
    def __init__(self, parent):
        super(CdhitFilterAgent, self).__init__(parent)
        options = [
            {"name": "query", "type": "infile", "format": "denovo_rna_v2.trinity_fasta"},  # 输入fasta文件
            {"name": "qunum", "type": "int", "default": 0},  # fasta编号
            {"name": "identity", "type": "float", "default": 0.99},  ##给出cdhit的参数identity
            {"name": "coverage", "type": "float", "default": 0.9},  # 给出cdhit的参数coverage
            {"name": "memory_limit", "type": "int", "default": 50000},  # 内存大小，0为无限制
            {"name": "method", "type": "int", "default": 0},  # 1为全局比对，0为局部比对
            {"name": "direction", "type": "int", "default": 1},  # 1为双向比对，0为单向比对
            {"name": "num_threads", "type": "int", "default": 10},  # cpu数
            {"name": "select", "type": "int", "default": 1},  # 1为聚类到最相似的类中，0为聚类到第一个符合阈值的类
            {"name": "out", "type": "string", "default": "cdhit.fasta"},  # 比对结果输出文件头
            {"name": "cdhit_filter_fasta", "type": "outfile", "format": "denovo_rna_v2.trinity_fasta"},  # 输入fasta文件
        ]
        self.add_option(options)
        self.step.add_steps('cdhitfilter')
        self.on('start', self.step_start)
        self.on('end', self.step_end)

    def step_start(self):
        self.step.cdhitfilter.start()
        self.step.update()

    def step_end(self):
        self.step.cdhitfilter.finish()
        self.step.update()

    def check_options(self):
        """
        重写参数检测函数
        """
        if not self.option("query").is_set:
            raise OptionError("必须设置参数query", code = "32002501")
        if not 0.75 <= self.option("identity") <= 1:
            raise OptionError("identity必须在0.75，1之间", code = "32002502")
        if not 0 <= self.option("coverage") <= 1:
            raise OptionError("coverage必须在0,1之间", code = "32002503")

    def set_resource(self):
        """
        设置所需资源
        """
        self._cpu = self.option("num_threads")
        self._memory = str(self.option("memory_limit") / 1000) + 'G'

    def end(self):
        result_dir = self.add_upload_dir(self.output_dir)
        result_dir.add_relpath_rules([
            [".", "", "结果输出目录"],
        ])
        super(CdhitFilterAgent, self).end()


class CdhitFilterTool(Tool):
    def __init__(self, config):
        super(CdhitFilterTool, self).__init__(config)
        self._version = '1.0'
        self.cdhit_est_path = 'bioinfo/uniGene/cd-hit-v4.5.7-2011-12-16/cd-hit-est'

    def run(self):
        super(CdhitFilterTool, self).run()
        self.run_cdhit()
        self.set_output()

    def word_len(self):
        word_length = 8
        if self.option("identity") >= 0.9:
            word_length = 8
        elif 0.88 <= self.option("identity") < 0.9:
            word_length = 7
        elif 0.85 <= self.option("identity") < 0.88:
            word_length = 6
        elif 0.8 <= self.option("identity") < 0.85:
            word_length = 5
        elif 0.75 <= self.option("identity") < 0.8:
            word_length = 4
        return word_length

    def run_cdhit(self):
        length = self.word_len()
        cmd = '%s -i %s -o %s -c %s -aS %s -n %s -G %s -M %s -d %s -r %s -g %s -T %s' % (
            self.cdhit_est_path, self.option("query").prop['path'], self.option("out"), self.option("identity"),
            self.option("coverage"), length, self.option("method"), self.option("memory_limit"), 0,
            self.option("direction"), self.option("select"), self.option("num_threads"))
        self.logger.info(cmd)
        command = self.add_command('cdhit', cmd)
        command.run()
        self.wait(command)
        if command.return_code == 0:
            self.logger.info("compare single succeed")
        else:
            self.set_error("compare single failed", code = "32002504")

    def set_output(self):
        os.link(self.option("out"), self.output_dir + "/cdhitfilter.fa")
        self.option("cdhit_filter_fasta", self.output_dir + "/cdhitfilter.fa")
        self.end()

    def linkdir(self, dirpath, dirname):
        """
        link一个文件夹下的所有文件到本module的output目录
        """
        allfiles = os.listdir(dirpath)
        newdir = os.path.join(self.output_dir, dirname)
        if not os.path.exists(newdir):
            os.mkdir(newdir)
        oldfiles = [os.path.join(dirpath, i) for i in allfiles]
        newfiles = [os.path.join(newdir, i) for i in allfiles]
        for newfile in newfiles:
            if os.path.exists(newfile):
                os.remove(newfile)
        for i in range(len(allfiles)):
            os.link(oldfiles[i], newfiles[i])

class TestFunction(unittest.TestCase):
    """
    This is test for the tool. Just run this script to do test.
    """
    def test(self):
        import random
        from mbio.workflows.single import SingleWorkflow
        from biocluster.wsheet import Sheet
        test_dir = '/mnt/ilustre/users/sanger-dev/sg-users/liubinxu/test_denovo'
        data = {
            "id": "Cdhitfilter" + str(random.randint(1, 10000)),
            "type": "tool",
            "name": "denovo_rna_v2.cdhit_filter",
            "instant": True,
            "options": dict(
                query=test_dir + "/" + "Trinity.fasta",
                num_threads="8",
                out="transcript"
            )
        }
        wsheet = Sheet(data=data)
        wf = SingleWorkflow(wsheet)
        wf.run()


if __name__ == '__main__':
    unittest.main()
