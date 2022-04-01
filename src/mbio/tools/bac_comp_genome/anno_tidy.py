# -*- coding: utf-8 -*-
# __author__ = 'hao.gao'
# last_modify: 2019.09.23

from biocluster.agent import Agent
from biocluster.tool import Tool
from biocluster.core.exceptions import OptionError
from mbio.packages.align.blast.xml2table import xml2table
import subprocess
import os


class AnnoTidyAgent(Agent):
    """
    比较基因组注释总览表，
    """
    def __init__(self, parent):
        super(AnnoTidyAgent, self).__init__(parent)
        options = [
            {"name": "sample", "type": "string"},  # 样品名
            {"name": "gff", "type": "infile", "format": "gene_structure.gff3"},  # 比基因组gff文件
            {"name": "cog", "type": "infile", "format": "sequence.profile_table"},  # 比基因组gff文件
            {"name": "kegg", "type": "infile", "format": "sequence.profile_table"},  # 比基因组gff文件
            {"name": "cazy", "type": "infile", "format": "sequence.profile_table"},  # 比基因组gff文件
            {"name": "card", "type": "infile", "format": "sequence.profile_table"},  # 比基因组gff文件
            {"name": "vfdb", "type": "infile", "format": "sequence.profile_table"},  # 比基因组gff文件
            {"name": "tcdb", "type": "infile", "format": "sequence.profile_table"},  # 比基因组gff文件
            {"name": "phi", "type": "infile", "format": "sequence.profile_table"},  # 比基因组gff文件
            {"name": "signalp", "type": "infile", "format": "sequence.profile_table"},  # 比基因组gff文件
            {"name": "signalp1", "type": "infile", "format": "sequence.profile_table"},  # 比基因组gff文件
            {"name": "tmhmm", "type": "infile", "format": "sequence.profile_table"},  # 比基因组gff文件
            {"name": "secretory", "type": "infile", "format": "sequence.profile_table"},  # 比基因组gff文件
            {"name": "out", "type": "outfile", "format": "sequence.profile_table"},  # 比基因组gff文件
        ]
        self.add_option(options)

    def check_options(self):
        for i in ["gff","cog","kegg",'cazy','card','vfdb','tcdb','phi','secretory','signalp','signalp1','tmhmm']:
            self.logger.info(i)
            if not self.option("{}".format(i)).is_set:
                raise OptionError("必须设置输入文件")
        return True

    def set_resource(self):
        self._cpu = 2
        self._memory = '10G'

    def end(self):
        result_dir = self.add_upload_dir(self.output_dir)
        result_dir.add_relpath_rules([
            [".", "", "结果输出目录"],
        ])
        super(AnnoTidyAgent, self).end()


class AnnoTidyTool(Tool):
    def __init__(self, config):
        super(AnnoTidyTool, self).__init__(config)
        self._version = "1.0"
        self.python_path = self.config.SOFTWARE_DIR + "/program/Python/bin/python"
        self.python_script = self.config.PACKAGE_DIR + '/bac_comp_genome/anno_tidy.py'
        self.sample =self.option('sample')
        self.gff = self.option('gff').prop['path']
        self.cog = self.option('cog').prop['path']
        self.kegg = self.option('kegg').prop['path']
        self.cazy = self.option('cazy').prop['path']
        self.vfdb = self.option('vfdb').prop['path']
        self.card = self.option('card').prop['path']
        self.phi = self.option('phi').prop['path']
        self.tcdb = self.option('tcdb').prop['path']
        self.signalp = self.option('signalp').prop['path']
        self.signalp1 = self.option('signalp1').prop['path']
        self.tmhmm = self.option('tmhmm').prop['path']
        self.secretory = self.option('secretory').prop['path']
        self.out = self.output_dir + "/"+ self.option("sample")

    def run(self):
        """
        运行
        :return:
        """
        super(AnnoTidyTool, self).run()
        self.run_anno()
        self.set_output()
        self.end()

    def run_anno(self):
        cmd = '{} {} -n {} -gff {} -cog {} -kegg {} -cazy {} -vfdb {} -card {} -phi {} -tcdb {} -signalp {} -signalps {} -tmhmm {} -secretory {} -o {}'.format(self.python_path, self.python_script,
                 self.sample, self.gff,self.cog,self.kegg,self.cazy,self.vfdb, self.card,self.phi,self.tcdb,self.signalp,self.signalp1,self.tmhmm,self.secretory,self.out)
        self.logger.info(cmd)
        try:
            subprocess.check_output(cmd, shell=True)
            self.logger.info('运行cog_anno完成')
        except subprocess.CalledProcessError:
            self.set_error('运行cog_anno出错')

    def set_output(self):
        if len(os.listdir(self.output_dir)) == 12:
            self.logger.info("output right")
        else:
            self.logger.info("output wrong")