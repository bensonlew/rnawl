# -*- coding: utf-8 -*-
# __author__ = 'gaohao'
# version 1.0
# last_modify: 20119.10.19

import os
import re, shutil
from biocluster.agent import Agent
from biocluster.tool import Tool
from biocluster.core.exceptions import OptionError
from mbio.packages.bac_comp_genome.common_function import link_file
import pandas as pd
import subprocess

class GenomeAaiAgent(Agent):
    """
    比较基因组的aai计算
    """
    def __init__(self, parent):
        super(GenomeAaiAgent, self).__init__(parent)
        options = [
            {"name": "seq_dir", "type": "infile", "format": "sequence.fasta_dir"},  # 输入参考序列文件夹
            {"name": "out", "type": "outfile", "format": "sequence.profile_table"},  #输出ani计算结果
            {"name": "evalue", "type": "string", "default": "1e-3"},
            {"name": "identity", "type": 'int', "default": 30},
            {"name": "aln_len", "type": 'int', "default": 70},
            {"name": "file_ext", "type": "string", "default": "fna"},
        ]
        self.add_option(options)

    def check_options(self):
        if not self.option("seq_dir").is_set:
            raise OptionError("必须设置参数seq_dir文件夹!")

    def set_resource(self):
        num = len(os.listdir(self.option("seq_dir").prop["path"]))
        self._cpu = 4
        self._memory = str(num*4)+'G'

    def end(self):
        super(GenomeAaiAgent, self).end()

class GenomeAaiTool(Tool):
    def __init__(self, config):
        super(GenomeAaiTool, self).__init__(config)
        self.ani = "program/Python/bin/comparem"
        path = os.path.join(self.config.SOFTWARE_DIR, "bioinfo/metaGenomic/Prodigal-2.6.3")
        self.set_environ(PATH=path)
        self.cmd_path = os.path.join(
            self.config.SOFTWARE_DIR, 'bioinfo/statistical/scripts/plot-hcluster_tree.pl')

    def run_aai(self):
        cmd = "{} aai_wf {} {} -e {} -p {} -a {} -x {} -c {}".format(self.ani, self.option("seq_dir").prop["path"], self.work_dir + "/out", self.option("evalue"), self.option("identity"), self.option("aln_len"), self.option("file_ext"), 4)
        command = self.add_command("run_aai", cmd).run()
        self.wait(command)
        if command.return_code == 0:
            self.logger.info("run_aai运行完成！")
        else:
            self.set_error("run_aai运行完成运行出错!")

    def set_output(self):
        link_file(self.work_dir + "/aai_summary.xls", self.output_dir + "/aai_summary.xls")
        self.option("out", self.output_dir + "/aai_summary.xls")

    def run(self):
        super(GenomeAaiTool, self).run()
        self.run_aai()
        self.summary_stat(self.work_dir + "/out/aai/aai_summary.tsv", self.work_dir + "/aai_summary.xls")
        self.set_output()
        self.end()

    def summary_stat(self, file, out):
        self.add_file(file,self.work_dir + "/all_result.xls")
        a = pd.read_table(self.work_dir + "/all_result.xls", sep='\t', header=0,)
        a = a.groupby(['#Genome A', 'Genome B'])['Mean AAI'].apply(lambda x:float(x)).unstack()
        a.index.name = "Genome A"
        a = a.fillna("100.0")
        a.to_csv(out, sep='\t', header=True, index=True)

    def add_file(self, file, file2):
        with open(file, "r") as f, open(file2, "w") as g:
            lines = f.readlines()
            g.write(lines[0])
            for line in lines[1:]:
                g.write(line)
                lin = line.strip().split("\t")
                g.write(lin[2] + "\t" + lin[1] + "\t" + lin[0] + "\t" + "\t".join(lin[3:]) + "\n")