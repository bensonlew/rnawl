#-*- coding: utf-8 -*-
# __author__ = 'qingchen.zhang'@20190326

import os
import re
from biocluster.agent import Agent
from biocluster.tool import Tool
from biocluster.core.exceptions import OptionError

class BamFastqAgent(Agent):
    """
    将上传的bam文件转为fastq文件
    """
    def __init__(self, parent):
        super(BamFastqAgent, self).__init__(parent)
        options = [
            {"name": "bam_file", "type": "infile", "format": "align.bwa.bam"},
        ]
        self.add_option(options)

    def check_options(self):
        if not self.option("bam_file").is_set:
            raise OptionError("必须添加bam_file文件！")

    def set_resource(self):
        self._cpu = 2
        self._memory = '10G'

    def end(self):
        super(BamFastqAgent, self).end()


class BamFastqTool(Tool):
    def __init__(self, config):
        super(BamFastqTool, self).__init__(config)
        self.samtools = "bioinfo/align/samtools-1.4/bin/samtools"

    def run_bam_fastq(self):
        """
        转格式bam>fastq
        :return:
        """
        bam_path = self.option('bam_file').prop['path']
        fastq1 = self.work_dir + '/all.l.fq'
        fastq2 = self.work_dir + '/all.r.fq'
        fastqs = self.work_dir + '/all.s.fq'
        cmd_fastq = "{} fastq -1 {} -2 {} -s {} --threads 4 {}".format(self.samtools, fastq1, fastq2, fastqs, bam_path)
        self.logger.info(cmd_fastq)
        to_fastq = 'to_fastq'
        command = self.add_command(to_fastq, cmd_fastq).run()
        self.wait(command)
        if command.return_code == 0:
            self.logger.info("%s运行完成"%to_fastq)
        else:
            self.set_error("%s运行失败"%to_fastq)

    def create_list(self):
        """
        生成list文件
        :return:
        """
        list_path = self.output_dir +"/list.txt"
        with open(list_path, 'w') as w:
            w.write('Sample\tSample_path\tType\n')
            w.write('all\tall.l.fq\tl\n')
            w.write('all\tall.r.fq\tr\n')
            if os.path.getsize(self.work_dir + '/all.s.fq') != 0:
                w.write('all\tall.s.fq\tr\n')

    def set_output(self):
        """
        设置结果目录
        :return:
        """
        if os.path.exists(self.output_dir + '/all.l.fq'):
            os.remove(self.output_dir + '/all.l.fq')
        os.link(self.work_dir + '/all.l.fq', self.output_dir + '/all.l.fq')
        if os.path.exists(self.output_dir + '/all.r.fq'):
            os.remove(self.output_dir + '/all.r.fq')
        os.link(self.work_dir + '/all.r.fq', self.output_dir + '/all.r.fq')
        if os.path.exists(self.output_dir + '/all.s.fq'):
            os.remove(self.output_dir + '/all.s.fq')
        if os.path.getsize(self.work_dir + '/all.s.fq') != 0:
            os.link(self.work_dir + '/all.s.fq', self.output_dir + '/all.s.fq')
        #self.create_list()

    def run(self):
        super(BamFastqTool, self).run()
        self.run_bam_fastq()
        self.set_output()
        self.end()