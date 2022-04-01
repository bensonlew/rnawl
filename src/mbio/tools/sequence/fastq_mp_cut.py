# -*- coding: utf-8 -*-
# __author__ = 'hao.gao'

import os, sys
# import shutil
from biocluster.core.exceptions import OptionError
from biocluster.agent import Agent
from biocluster.tool import Tool


class FastqMpCutAgent(Agent):
    """
    细菌基因组MP文库fq序列截取
    version: v1
    author: gaohao
    last_modify: 2018.01.25
    """

    def __init__(self, parent):
        super(FastqMpCutAgent, self).__init__(parent)
        options = [
            {"name": "fastq1", "type": "infile", "format": "sequence.fastq"},  # 输入文件,sample.mp.sickle.l.fastq
            {"name": "fastq2", "type": "infile", "format": "sequence.fastq"},  # 输入文件,sample.mp.sickle.r.fastq
            {"name": "cut_num", "type": "int", "default": 40},  # 截取长度，默认40bp
            {'name': 'sample_name', "type": "string"},  # 样本名
        ]
        self.add_option(options)


    def check_options(self):
        """
        检查参数
        :return:
        """
        if not self.option('fastq1'):
            raise OptionError('必须输入*l.fastq文件', code="34001001")
        if not self.option('fastq2'):
            raise OptionError('必须输入*2.fastq文件', code="34001002")

    def set_resource(self):
        """
        :return:
        """
        self._cpu = 2
        self._memory = "5G"

    def end(self):
        super(FastqMpCutAgent, self).end()


class FastqMpCutTool(Tool):
    def __init__(self, config):
        super(FastqMpCutTool, self).__init__(config)
        self.fq1 = self.option('fastq1').prop['path']
        self.fq2 = self.option('fastq2').prop['path']
        self.num = self.option('cut_num')
        self.sample_name =self.option('sample_name')
        self.perl_path = "/program/perl/perls/perl-5.24.0/bin/perl"
        self.perl_script = self.config.PACKAGE_DIR + "/bacgenome/"


    def run(self):
        """
        运行
        :return:
        """
        super(FastqMpCutTool, self).run()
        self.cut_fq()
        self.set_output()
        self.end()

    def cut_fq(self):
        num =0
        for file in [self.fq1,self.fq2]:
            num += 1
            outfile = self.work_dir + '/' + self.sample_name + '.trim-cut' + str(self.num) + '.' + str(num) + '.fq'
            cmd = '{} {}fastq_cut.pl {} {} {}'.format(self.perl_path,self.perl_script,file,self.num,outfile)
            run_fa = 'cut_fq' + str(num)
            command = self.add_command(run_fa, cmd).run()
            self.wait(command)
            if command.return_code == 0:
                self.logger.info("%s运行完成" % run_fa)
            else:
                self.set_error("%s运行出错!" , variables=( run_fa), code="34001001")

    def set_output(self):
        """
        将结果文件复制到output文件夹下面
        :return:
        """
        self.logger.info("设置MP序列截取目录")
        file1 = self.output_dir + '/' + self.sample_name + '.trim-cut' + str(self.num) + '.' + '1.fq'
        file2 = self.output_dir + '/' + self.sample_name + '.trim-cut' + str(self.num) + '.' + '2.fq'
        for file in[file1,file2]:
            if os.path.exists(file):
                os.remove(file)
        os.link(self.work_dir + '/' + self.sample_name + '.trim-cut' + str(self.num) + '.' + '1.fq',file1)
        os.link(self.work_dir + '/' + self.sample_name + '.trim-cut' + str(self.num) + '.' + '2.fq', file2)
        self.logger.info("设置MP截取序列结果目录成功")
