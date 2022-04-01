# -*- coding: utf-8 -*-
#__author__ = 'hao.gao'

import os
import re
from biocluster.agent import Agent
from biocluster.tool import Tool
from biocluster.core.exceptions import OptionError
from mbio.packages.metagbin.common_function import link_file

class InquireSeqAgent(Agent):
    """
    序列查询的blastn or blastx
    """
    def __init__(self, parent):
        super(InquireSeqAgent, self).__init__(parent)
        options = [
            {"name": "in_file", "type": "infile", "format": "sequence.fasta"}, #
            {"name": "method", "type": "string"},  #主要是blastn or blastx
            {"name": "ref", "type": "infile", "format": "sequence.fasta"},  #参考数据库
        ]
        self.add_option(options)

    def check_options(self):
        if not self.option("in_file").is_set:
            raise OptionError("必须添加查询序列文件！")
        if not self.option("method"):
            raise OptionError("必须添加查询方法！")
        else:
            if self.option("method") not in ['blastn','blastx']:
                raise OptionError("必须添加正确的查询方法！")
        if not self.option("ref").is_set:
            raise OptionError("必须添加参考序列文件！")

    def set_resource(self):
        self._cpu = 2
        self._memory = '10G'

    def end(self):
        super(InquireSeqAgent, self).end()

class InquireSeqTool(Tool):
    def __init__(self, config):
        super(InquireSeqTool, self).__init__(config)
        self.ref = self.option("ref").prop['path']
        self.blast = "/bioinfo/align/ncbi-blast-2.3.0+/bin/"
        self.out = self.work_dir + '/all.blast.xls'

    def run_makedb(self):
        """
        将sort结果进行merge
        :return:
        """
        dbtype = ''
        if self.option("method") in ['blastn']:
            dbtype='nucl'
        elif self.option("method") in ['blastx']:
            dbtype = 'prot'
        cmd = "{}makeblastdb -in {} -dbtype {} -out {}".format(self.blast, self.ref, dbtype,self.work_dir + "/db")
        self.logger.info(cmd)
        command1 = self.add_command('run_makedb', cmd).run()
        self.wait(command1)
        if command1.return_code == 0:
            self.logger.info("run_makedb运行完成")
        else:
            self.set_error("run_makedb运行失败")

    def run_blast(self):
        """
        序列比对
        :return:
        """
        self.input = self.option("in_file").prop['path']
        type = ''
        if self.option("method") in ['blastn']:
            type='blastn'
        elif self.option("method") in ['blastx']:
            type = 'blastx'
        cmd = "{}{} -query {} -db {} -out {} -evalue 1e-5 -max_target_seqs 1 -num_threads 4 -outfmt 6".format(self.blast,type,self.input,self.work_dir + '/db',self.out)
        self.logger.info(cmd)
        command = self.add_command('run_blast', cmd).run()
        self.wait(command)
        if command.return_code == 0:
            self.logger.info("run_blast运行完成")
        else:
            self.set_error("run_blast运行失败")

    def set_output(self):
        """
        设置结果文件目录
        :return:
        """
        link_file(self.out,self.output_dir + '/all.blast.xls')

    def run(self):
        super(InquireSeqTool, self).run()
        self.run_makedb()
        self.run_blast()
        self.set_output()
        self.end()