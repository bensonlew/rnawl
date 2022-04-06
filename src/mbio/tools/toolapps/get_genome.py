# -*- coding: utf-8 -*-
# __author__ = 'hao.gao'

from biocluster.core.exceptions import OptionError
from biocluster.agent import Agent
from biocluster.tool import Tool
import pandas as pd
import os

class GetGenomeAgent(Agent):
    '''
    16s比对数据库（GTDB），提取数据
    '''
    def __init__(self, parent):
        super(GetGenomeAgent, self).__init__(parent)
        options = [
            {'name': 'blast', 'type': 'infile', "format": "sequence.profile_table"},  ## 16s比对的结果
            {"name": "ref_s16", "type": "outfile", "format": "sequence.fasta"},
            {"name": "ref_genome", "type": "outfile", "format": "sequence.fasta_dir"},
            {"name": "cus_table", "type": "infile", "format": "sequence.profile_table"},
            {'name': 'blast_out', 'type': 'outfile', "format": "sequence.profile_table"}, ##按照98.7的indentity判断获取
            {"name": "s16_fa", "type": "outfile", "format": "sequence.fasta"},
            {"name": "genomes", "type": "outfile", "format": "sequence.fasta_dir"}
        ]
        self.add_option(options)

    def check_options(self):
        """
        参数二次检查
        :return:
        """
        if not self.option('blast').is_set:
            raise OptionError("必须输入blast!")
        return True

    def set_resource(self):
        """
        设置资源
        :return:
        """
        self._cpu = '2'
        self._memory = '10G'

    def end(self):
        """
        运行结束
        :return:
        """
        super(GetGenomeAgent, self).end()

class GetGenomeTool(Tool):
    """
    16s比对数据库（GTDB）
    """
    def __init__(self, config):
        super(GetGenomeTool, self).__init__(config)
        self.python = "/miniconda2/bin/python"
        self.script = self.config.PACKAGE_DIR + "/toolapps/"
        self.gtdb_16s = self.config.SOFTWARE_DIR + "/database/GTDB/GTDB_rrna/bac120_ssu_reps_r95.fna"
        self.gtdb_genomes = self.config.SOFTWARE_DIR + "/database/GTDB/release95/fastani/database/"


    def run(self):
        """
        运行
        :return:
        """
        super(GetGenomeTool, self).run()
        if self.option("cus_table").is_set:
            self.run_getcustomfa()
            self.set_output()
        else:
            self.run_getfasta()
            self.set_output()
        self.end()

    def run_getfasta(self):
        """
        获取数据
        :return:
        """
        cmd = '{} {}get_16s_genome.py -b {} -r {} -d {} -o {}'.format(self.python, self.script, self.option("blast").prop['path'], self.gtdb_16s, self.gtdb_genomes, self.output_dir)
        self.logger.info(cmd)
        command = self.add_command("run_getfasta", cmd)
        command.run()
        self.wait()
        if command.return_code == 0:
            self.logger.info('run_getfasta运行成功！')
        else:
            self.set_error('run_getfasta运行失败！')

    def run_getcustomfa(self):
        """
        获取数据
        :return:
        """
        cmd = '{} {}get_16s_custom.py -b {} -r {} -o {}'.format(self.python, self.script, self.option("blast").prop['path'], self.option("cus_table").prop['path'], self.output_dir)
        self.logger.info(cmd)
        command = self.add_command("run_getfasta", cmd)
        command.run()
        self.wait()
        if command.return_code == 0:
            self.logger.info('run_getfasta运行成功！')
        else:
            self.set_error('run_getfasta运行失败！')

    def set_output(self):
        """
        设置结果目录和结果文件
        :return:
        """
        if os.path.exists(self.output_dir+"/all.blast_ani.xls"):
            self.option("blast_out", self.output_dir+"/all.blast_ani.xls")
        if os.path.getsize(self.output_dir+"/all.16s.fa") >0:
            self.option("s16_fa", self.output_dir+"/all.16s.fa")
        if os.path.getsize(self.output_dir+"/genome") >0:
            self.option("genomes", self.output_dir+"/genome")