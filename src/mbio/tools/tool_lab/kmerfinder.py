# -*- coding: utf-8 -*-
# __author__ = 'gaohao'
# version 1.0
# last_modify: 2021.04.29

import os, time
import re, shutil
from biocluster.agent import Agent
from biocluster.tool import Tool
from biocluster.core.exceptions import OptionError
from mbio.packages.bac_comp_genome.common_function import link_dir
import json
import pandas as pd

class KmerfinderAgent(Agent):
    """
    单个基因组Kmerfinder预测
    """
    def __init__(self, parent):
        super(KmerfinderAgent, self).__init__(parent)
        options = [
            {"name": "fasta", "type": "infile", "format": "sequence.fasta"},# 序列文件
            {"name": "sample_name", "type": "string"},# 样本名
            {"name": "species", "type": "string", },# 物种名称
            {"name": "method", "type": "string", 'default':"blastn"},  #方法 kmer or blastn
        ]
        self.add_option(options)

    def check_options(self):
        if not self.option("fasta").is_set:
            raise OptionError("必须输入序列文件")
        if self.option('species') not in ['viral', 'bacteria', 'fungi', 'typestrain', 'archaea', 'protozoa']:
            raise OptionError("必须输入正确的来源信息！")

    def set_resource(self):
        self._cpu = 2
        self._memory = '20G'

    def end(self):
        super(KmerfinderAgent, self).end()

class KmerfinderTool(Tool):
    def __init__(self, config):
        super(KmerfinderTool, self).__init__(config)
        self.path = self.config.SOFTWARE_DIR + "/program/Python35/bin"
        self.lib = self.config.SOFTWARE_DIR + "/program/Python35/lib"
        self.set_environ(PATH=self.path, LD_LIBRARY_PATH=self.lib)
        self.python = "/miniconda2/bin/python"
        self.python2 = "/program/Python35/bin/python"
        self.python_script = self.config.SOFTWARE_DIR + "/bioinfo/Genomic/Sofware/kmerfinder_v3/"
        self.blastn = self.config.SOFTWARE_DIR + "/bioinfo/align/ncbi-blast-2.3.0+/bin/blastn"
        self.kma = self.config.SOFTWARE_DIR + "/bioinfo/Genomic/Sofware/kma/"
        if self.option('species') in ['viral']:
            self.db = self.config.SOFTWARE_DIR + "/bioinfo/Genomic/Sofware/kmerfinder_v3/kmerfinder_db/viral/viral.ATG"
            self.tax = self.config.SOFTWARE_DIR + "/bioinfo/Genomic/Sofware/kmerfinder_v3/kmerfinder_db/viral/viral.tax"
        elif self.option('species') in ['bacteria']:
            self.db = self.config.SOFTWARE_DIR + "/bioinfo/Genomic/Sofware/kmerfinder_v3/kmerfinder_db/bacteria/bacteria.ATG"
            self.tax = self.config.SOFTWARE_DIR + "/bioinfo/Genomic/Sofware/kmerfinder_v3/kmerfinder_db/bacteria/bacteria.tax"
        elif self.option('species') in ['fungi']:
            self.db = self.config.SOFTWARE_DIR + "/bioinfo/Genomic/Sofware/kmerfinder_v3/kmerfinder_db/fungi/fungi.ATG"
            self.tax = self.config.SOFTWARE_DIR + "/bioinfo/Genomic/Sofware/kmerfinder_v3/kmerfinder_db/fungi/fungi.tax"
        elif self.option('species') in ['typestrain']:
            self.db = self.config.SOFTWARE_DIR + "/bioinfo/Genomic/Sofware/kmerfinder_v3/kmerfinder_db/typestrain/typestrain.ATG"
            self.tax = self.config.SOFTWARE_DIR + "/bioinfo/Genomic/Sofware/kmerfinder_v3/kmerfinder_db/typestrain/typestrain.tax"
        elif self.option('species') in ['archaea']:
            self.db = self.config.SOFTWARE_DIR + "/bioinfo/Genomic/Sofware/kmerfinder_v3/kmerfinder_db/archaea/archaea.ATG"
            self.tax = self.config.SOFTWARE_DIR + "/bioinfo/Genomic/Sofware/kmerfinder_v3/kmerfinder_db/archaea/archaea.tax"
        elif self.option('species') in ['protozoa']:
            self.db = self.config.SOFTWARE_DIR + "/bioinfo/Genomic/Sofware/kmerfinder_v3/kmerfinder_db/protozoa/protozoa.ATG"
            self.tax = self.config.SOFTWARE_DIR + "/bioinfo/Genomic/Sofware/kmerfinder_v3/kmerfinder_db/protozoa/protozoa.tax"

    def run_kmerfinder(self):
        """
        运行mlst
        :return:
        """
        if os.path.exists(self.work_dir+"/result"):
           shutil.rmtree(self.work_dir+"/result")
        os.mkdir(self.work_dir+"/result")
        cmd = "{} {}kmerfinder.py -i {} -o {} -tax {} -kp {} -db {}".format(self.python2, self.python_script, self.option("fasta").prop['path'], self.work_dir+"/result", self.tax, self.kma, self.db)
        self.logger.info(cmd)
        command = self.add_command("run_kmerfinder", cmd).run()
        self.wait(command)
        if command.return_code == 0:
            self.logger.info("run_kmerfinder运行完成！")
        else:
            self.set_error("run_kmerfinder运行出错!")

    def run_stat(self):
        """
        运行处理文件，生成结果文件
        :return:
        """
        table = pd.read_table(self.work_dir+"/result/results.txt", sep="\t", header=0)
        table = table[['# Assembly', 'Score', 'tot_depth', 'q_value', 'p_value', 'Description', 'Accession Number', 'TAXID', 'Taxonomy']]
        table['Sample'] = self.option("sample_name")
        table.to_csv(self.output_dir + "/"+self.option("sample_name")+".KmerFinderResult.xls",sep="\t", header=True, index=False,columns=['Sample', '# Assembly', 'Score', 'tot_depth', 'q_value', 'p_value', 'Description', 'Accession Number', 'TAXID', 'Taxonomy'])


    def set_output(self):
        """
        设置结果文件目录
        """
        self.logger.info("生成结果文件完成")

    def run(self):
        """
        运行
        """
        super(KmerfinderTool, self).run()
        self.run_kmerfinder()
        if os.path.getsize(self.work_dir+"/result/results.txt") >0:
            self.run_stat()
        else:
            self.end()
        self.set_output()
        self.end()