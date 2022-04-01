# -*- coding: utf-8 -*-
# __author__ = 'zhangyitong'
import os,glob
import shutil
from biocluster.core.exceptions import OptionError
from biocluster.agent import Agent
from biocluster.tool import Tool
from biocluster.config import Config
from mbio.packages.whole_transcriptome.utils import runcmd
import unittest
import json
from collections import OrderedDict
from Bio import SeqIO
import pandas as pd


class ToolGenesetPpiAgent(Agent):
    """
    This script is used to extracts mapped/unmapped reads from input bam/sam file.
    """
    def __init__(self, parent):
        super(ToolGenesetPpiAgent, self).__init__(parent)
        options = [
            {"name": "seq_file", "type": "infile", "format": "ref_rna_v2.fasta"},
            {"name": "geneset", "type": "infile", "format": 'ref_rna_v2.common'},
        ]
        self.add_option(options)
        self.step.add_steps("extract")
        self.on('start', self.stepstart)
        self.on('end', self.stepfinish)

    def stepstart(self):
        self.step.extract.start()
        self.step.update()

    def stepfinish(self):
        self.step.extract.finish()
        self.step.update()

    def check_options(self):
        return True

    def set_resource(self):
        """
        设置所需资源，需在之类中重写此方法 self._cpu ,self._memory
        :return:
        """
        self._cpu = 1
        self._memory = "20G"

    def end(self):
        super(ToolGenesetPpiAgent, self).end()


class ToolGenesetPpiTool(Tool):
    def __init__(self, config):
        super(ToolGenesetPpiTool, self).__init__(config)
        self.set_environ(PATH=os.path.join(self.config.SOFTWARE_DIR, 'program/Python/bin'))

    def run(self):
        """
        运行
        :return:
        """
        super(ToolGenesetPpiTool, self).run()
        self.run_tool()
        self.set_output()
        self.end()

    def run_tool(self):
        cds_fa_dict = dict()
        record_fa = SeqIO.parse(self.option('seq_file').prop['path'], 'fasta')
        for i in record_fa:
            cds_fa_dict[i.id] = str(i.seq)

        geneset_df = pd.read_table(self.option('geneset').prop['path'], header=0)
        geneset_df['seq'] = geneset_df.apply(lambda x: cds_fa_dict.get(x['gene_id']), axis=1)
        geneset_df.to_csv(os.path.join(self.work_dir, 'geneset_file.txt'), sep='\t', index=False)

    def set_output(self):
        """
        将结果文件复制到output文件夹下面
        :return:
        """
        self.logger.info("设置结果目录")
