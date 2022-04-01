# -*- coding: utf-8 -*-
# __author__ = 'qindanhua'
from biocluster.agent import Agent
from biocluster.tool import Tool
import os
from biocluster.core.exceptions import OptionError
import subprocess


class FastaStatAgent(Agent):
    """
    fastastat:用于统计fasta文件中序列数目，长度等信息
    version 1.0
    author: qindanhua
    last_modify: 2016.01.05
    """

    def __init__(self, parent):
        super(FastaStatAgent, self).__init__(parent)
        options = [
            {"name": "fasta", "type": "infile", "format": "sequence.fasta"}
        ]
        self.add_option(options)
        self.step.add_steps('seqstat')
        self.on('start', self.step_start)
        self.on('end', self.step_end)

    def step_start(self):
        self.step.seqstat.start()
        self.step.update()

    def step_end(self):
        self.step.seqstat.finish()
        self.step.update()

    def check_options(self):
        """
        检查参数是否正确
        """
        if not self.option("fasta").is_set:
            raise OptionError("请传入OTU代表序列文件")

    def set_resource(self):
        """
        所需资源
        """
        self._cpu = 10
        self._memory = ''


class FastaStatTool(Tool):
    """
    version 1.0
    """

    def __init__(self, config):
        super(FastaStatTool, self).__init__(config)
        self.seqstat_path = 'bioinfo/seq/scripts/'
        # self.python_path = '/program/Python/bin/'

    def seqstat(self):
        """
        运行seqstat工具，统计fasta文件的序列长度等信息
        :return:
        """
        cmd = self.seqstat_path + ' -i %s' % (self.seqstat_path, self.option('fasta').prop['path'])
        self.logger.info('开始统计fasta文件信息')
        print self.config.SOFTWARE_DIR + cmd
        try:
            subprocess.check_output(self.config.SOFTWARE_DIR + cmd, shell=True)
            self.logger.info("统计完成")
            return True
        except subprocess.CalledProcessError:
            self.logger.info("统计过程出现错误")
            return False

    def set_output(self):
        """
        设置输出文件路径
        :return:
        """
        self.logger.info("set output")
        for f in os.listdir(self.output_dir):
            os.remove(os.path.join(self.output_dir, f))
        os.link(self.work_dir+'/fasta_stat.xls', self.output_dir+'/fasta_stat.xls')
        os.link(self.work_dir+'/fasta_len.xls', self.output_dir+'/fasta_len.xls')
        self.logger.info("done")

    def run(self):
        super(FastaStatTool, self).run()
        self.seqstat()
        self.set_output()
        self.end()
