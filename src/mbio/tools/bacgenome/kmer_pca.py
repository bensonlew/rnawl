# -*- coding: utf-8 -*-
# __author__ = 'haidong.gu'
# __last_modify__ = '2019/4/16'
from biocluster.agent import Agent
from biocluster.tool import Tool
import os,sys
from biocluster.core.exceptions import OptionError
from mbio.packages.metagenomic.common import link_dir
import pandas as pd


class KmerPcaAgent(Agent):
    """
    DemoAgent:
    version 1.0
    """

    def __init__(self, parent):
        super(KmerPcaAgent, self).__init__(parent)
        options = [
            {"name": "scaf", "type": "infile", "format": "sequence.fasta", "required": True},
            {"name": "min_len", "type": "int", "default": 1000, "required": True},
            {"name": "kmer", "type": "int", "default": 4, "required": True},
            {"name": "cov", "type": "infile", "format": "sequence.profile_table"}  # coverage信息
        ]
        self.add_option(options)

    def check_options(self):
        """
        检查参数是否正确
        """
        # modified check_option

    def set_resource(self):
        """
        所需资源
        """
        self._cpu = 2
        self._memory = "10G"


class KmerPcaTool(Tool):
    def __init__(self, config):
        super(KmerPcaTool, self).__init__(config)
        self.python_path = self.config.SOFTWARE_DIR + '/miniconda2/bin/python'
        self.perl_path = "/program/perl-5.24.0/bin/perl"
        self.R_path = 'program/R-3.3.1/bin/Rscript'
        self.script_path = self.config.SOFTWARE_DIR + "/bioinfo/metaGenomic/CONCOCT-develop/scripts/cut_up_fasta.py"
        self.perl_script = self.config.PACKAGE_DIR + "/statistical/ordination.pl"
        kmer_cal_packages = self.config.SOFTWARE_DIR + "/bioinfo/metaGenomic/CONCOCT-develop/concoct"
        sys.path.append(kmer_cal_packages)
        self.split_fasta = os.path.join(self.work_dir, "split.fasta")
        # self.freq_table = os.path.join(self.work_dir, "output.xls")
        self.freq_table = os.path.join(self.work_dir, "%smer" % self.option("kmer"))

    def run_split(self):
        """
        description
        :return:
        """
        cmd = "%s %s -c %s -o 0 %s > %s" % (self.python_path, self.script_path, self.option("min_len"), self.option("scaf").prop["path"], self.split_fasta)
        command = self.add_command('split_fasta', cmd, shell=True).run()
        self.wait(command)
        if command.return_code == 0:
            self.logger.info("截取序列成功")
        else:
            self.set_error("截取序列失败")

    def run_calculate(self):
        from input import load_composition
        comp, length = load_composition(self.split_fasta, self.option("kmer"), self.option("min_len")-1)
        comp.T.to_csv(self.freq_table, sep="\t", index=True)

    def run_pca(self):
        # perl ~/biocluster/src/mbio/packages/statistical/ordination.pl -type pca -community output.xls -outdir . -scale F
        cmd = "%s %s -type pca -community %s -outdir . -scale F" % (self.perl_path, self.perl_script, self.freq_table)
        command = self.add_command("pca", cmd).run()
        self.wait(command)
        if command.return_code == 0:
            self.logger.info("pca成功")
        else:
            self.set_error("pca失败")

    def run_r(self):
        cmd = self.R_path + " cmd.r"
        command = self.add_command("pca_r", cmd).run()
        self.wait(command)
        if command.return_code == 0:
            self.logger.info("r成功")
        else:
            self.set_error("r失败")

    def get_cov(self, df, cov):
        seq_id = df.name
        seq_id = os.path.splitext(seq_id)[0]
        cov_value = cov.loc[seq_id,"totalAvgDepth"]
        return cov_value

    def set_output(self):
        """
        设置输出文件路径
        :return:
        """
        if self.option("cov").is_set:
            self.logger.info("have coverage")
            cov_map = pd.read_table(self.option("cov").prop["path"], index_col="contigName")
            pca_data = pd.read_table(self.work_dir + "/pca/" + os.path.basename(self.freq_table) + "_pca_sites.xls",index_col="Sample_ID")
            pca_data["cov"] = pca_data.apply(self.get_cov, axis=1, cov=cov_map)
            pca_data.reindex(columns=["PC1", "PC2", "cov"]).to_csv(self.work_dir + "/pca/" + os.path.basename(self.freq_table) + "_pca_cov_sites.xls", sep="\t", header=True, index=True)
        else:
            self.logger.info("not have coverage")
        link_dir(os.path.join(self.work_dir, "pca"), self.output_dir)


    def run(self):
        super(KmerPcaTool, self).run()
        self.run_split()
        self.run_calculate()
        self.run_pca()
        self.run_r()
        self.set_output()
        self.end()