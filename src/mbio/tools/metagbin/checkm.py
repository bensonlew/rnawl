#!/usr/bin/env python
# -*- coding: utf-8 -*-
# __author__ = 'gaohao'

from biocluster.agent import Agent
from biocluster.tool import Tool
from biocluster.core.exceptions import OptionError
import os,shutil

class CheckmAgent(Agent):
    """
    用于评估bin和genome的完整性
    version 1.0
    author: gaohao
    last_modify: 2019.01.07
    """

    def __init__(self, parent):
        super(CheckmAgent, self).__init__(parent)
        options = [
            {"name": "bin_dir", "type": "infile", "format": "sequence.fasta_dir"},  #bin的文件目录
            {"name": "method", "type": "string", "default": "bin_checkm"},  ##组装序列
        ]
        self.add_option(options)
        self.list =[]


    def check_options(self):
        """
        检测参数是否正确
        """
        if not self.option('bin_dir').is_set:
            raise OptionError("生成bin的fasta文件夹不存在！")


    def set_resource(self):
        """
        所需资源
        """
        files =os.listdir(self.option('bin_dir').prop['path'])
        if len(files) ==1:
            num =int(len(files)*20)
            self._cpu = 4
            self._memory = str(num) + 'G'
        elif len(files) > 1:
            num = int(len(files) * 2.5)
            self._cpu = 8
            self._memory = str(num) + 'G'

    def end(self):
        super(CheckmAgent, self).end()


class CheckmTool(Tool):
    """
    version 1.0
    """
    def __init__(self, config):
        super(CheckmTool, self).__init__(config)
        self.path =self.config.SOFTWARE_DIR + "/bioinfo/metaGenomic/CheckM-master:" + self.config.SOFTWARE_DIR + "/miniconda2/bin:" + self.config.SOFTWARE_DIR + "/bioinfo/metaGenomic/Prodigal-2.6.3:" + self.config.SOFTWARE_DIR + "/bioinfo/Genomic/Sofware/islandpath_dimob/hmmer-3.1b1/bin:" + self.config.SOFTWARE_DIR + "/bioinfo/metaGenomic/pplacer-Linux-v1.1.alpha19"
        self.set_environ(PATH=self.path)
        self.checkm = "/miniconda2/bin/"
        self.perl = "/miniconda2/bin/perl"
        self.perl_script = self.config.PACKAGE_DIR + "/metagbin/checkm_stat.pl"
        self.bin_dir =self.option('bin_dir').prop['path']
        self.summary =self.work_dir + '/checkm.summary.txt'
        self.result = self.work_dir + '/result'
        self.tetra_result = self.work_dir + '/tetra.summary.xls'
        self.taxon = self.config.SOFTWARE_DIR + "/bioinfo/metaGenomic/CheckM-master/checkm.taxon.xls"

    def run_checkm(self):
        if os.path.exists(self.result):
            shutil.rmtree(self.result)
        cmd = "{}checkm lineage_wf -x fa -t 16 --pplacer_threads 16 -f {} --tab_table --tmpdir {} {} {}".format(self.checkm,self.summary, "./",self.bin_dir,self.result)
        self.logger.info(cmd)
        self.logger.info("开始运行run_checkm")
        command = self.add_command("run_checkm", cmd)
        command.run()
        self.wait(command)
        if command.return_code == 0:
            self.logger.info("运行run_checkm完成")
        else:
            self.set_error("运行run_checkm运行出错!")

    def run_checkm_genome(self):
        """
        古菌基因组评估
        """
        self.summary =self.work_dir + '/checkm.summary.txt'
        out_file = self.output_dir + '/checkm_assess.xls'
        with open(self.summary, 'r') as infile, open(out_file, 'w') as outfile:
            outfile.write('#Completeness\tContamination\tStrain heterogeneity\n')
            for line in infile[1:]:
                line = line.strip().split('\t')
                complete = line[-3]
                contamina = line[-2]
                heter = line[-1]
                outfile.write('{}\t{}\t{}\n'.format(complete, contamina, heter))
        self.logger.info("古菌评估文件写入成功")

    def run_checkm_stat(self):
        cmd = "{} {} {} {} {} {}".format(self.perl,self.perl_script,self.summary,self.result + '/storage/bin_stats_ext.tsv',self.taxon,"all")
        self.logger.info(cmd)
        self.logger.info("开始运行run_checkm_stat")
        command = self.add_command("run_checkm_stat", cmd)
        command.run()
        self.wait(command)
        if command.return_code == 0:
            self.logger.info("运行run_checkm_stat完成")
        else:
            self.set_error("运行run_checkm_stat运行出错!")

    def set_output(self):
        for i in ['all.marker.xls','all.bin.summary.xls']:
            if os.path.exists(self.output_dir + "/" + i):
                os.remove(self.output_dir + "/" + i)
            os.link(self.work_dir + "/" + i,self.output_dir + "/" + i)

    def run(self):
        """
        运行
        """
        super(CheckmTool, self).run()
        if self.option("method") in ["bin_checkm"]:
            self.run_checkm()
            self.run_checkm_stat()
            self.set_output()
            self.end()
        else:
            self.run_checkm()
            self.run_checkm_genome()
            self.end()