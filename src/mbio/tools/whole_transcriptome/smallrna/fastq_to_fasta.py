# -*- coding: utf-8 -*-
# __author__ = 'qiuping'
from biocluster.agent import Agent
from biocluster.tool import Tool
from biocluster.core.exceptions import OptionError
import os

class FastqToFastaAgent(Agent):
    """
    将fastq文件转换成fasta格式的文件
    version v1
    author：qiuping
    last_modify:2015.01.06
    """
    def __init__(self, parent):
        super(FastqToFastaAgent, self).__init__(parent)
        options = [
            {"name": "fastq_input", "type": "infile", "format": "sequence.fastq"},
            {"name": "fasta_id", "type": "string", "default": "none"}
        ]
        self.add_option(options)
        self.step.add_steps("fastq_to_fasta")
        self.on('start', self.stepstart)
        self.on('end', self.stepfinish)

    def stepstart(self):
        self.step.fastq_to_fasta.start()
        self.step.update()

    def stepfinish(self):
        self.step.fastq_to_fasta.finish()
        self.step.update()

    def check_options(self):
        """
        检查参数设置
        :return:
        """
        if not self.option('fastq_input'):
            raise OptionError("必须设置输入的fastq文件")

    def set_resource(self):
        """
        设置所需资源
        :return:
        """
        self._cpu = 10
        self._memory = '5G'


class FastqToFastaTool(Tool):
    def __init__(self, config):
        super(FastqToFastaTool, self).__init__(config)
        self._version = "v1"
        self.cmd = "/miniconda2/bin/"
        self.script = self.config.SOFTWARE_DIR + "/bioinfo/seq/scripts/"

    def run(self):
        """
        运行
        :return:
        """
        super(FastqToFastaTool, self).run()
        self.fastq_to_fasta()
        self.end()

    def fastq_to_fasta(self):
        """
        运行fastq_to_fasta.py程序，将单个fastq文件转换成fasta文件
        """
        self.logger.info("开始运行pair_fastq_to_fasta命令")
        cmd = self.cmd + "python %sfastq_to_fasta.py -i %s -o fasta -n %s" % \
                         (self.script, self.option('fastq_input').prop['path'], self.option('fasta_id'))
        self.logger.info(cmd)
        fasta_command = self.add_command("fastq_to_fasta_cmd", cmd).run()
        self.wait(fasta_command)
        if fasta_command.return_code == 0:
            self.logger.info("fastq_to_fasta_cmd运行完成")
        else:
            self.set_error("fastq_to_fasta_cmd运行出错!")
        for root, dirs, files in os.walk(self.output_dir):
            for names in files:
                os.remove(os.path.join(root, names))
        os.link(self.work_dir + '/fasta', self.output_dir + '/fasta')
