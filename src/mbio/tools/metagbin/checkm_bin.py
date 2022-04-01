#!/usr/bin/env python
# -*- coding: utf-8 -*-
# __author__ = 'gaohao'

from biocluster.agent import Agent
from biocluster.tool import Tool
from biocluster.core.exceptions import OptionError
import os,shutil

class CheckmBinAgent(Agent):
    """
    用于评估bin
    version 1.0
    author: gaohao
    last_modify: 2019.03.07
    """

    def __init__(self, parent):
        super(CheckmBinAgent, self).__init__(parent)
        options = [
            {"name": "bin_dir", "type": "infile", "format": "sequence.fasta_dir"},  #bin的文件目录
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
            num =int(len(files)*50)
            self._cpu = 4
            self._memory = str(num) + 'G'
        elif len(files) > 1:
            num = int(len(files) * 2.5)
            self._cpu = 8
            self._memory = str(num) + 'G'

    def end(self):
        super(CheckmBinAgent, self).end()

class CheckmBinTool(Tool):
    """
    version 1.0
    """
    def __init__(self, config):
        super(CheckmBinTool, self).__init__(config)
        self.path =self.config.SOFTWARE_DIR + "/bioinfo/metaGenomic/CheckM-master:" + self.config.SOFTWARE_DIR + "/program/Python/bin:" + self.config.SOFTWARE_DIR + "/bioinfo/metaGenomic/Prodigal-2.6.3:" + self.config.SOFTWARE_DIR + "/bioinfo/Genomic/Sofware/islandpath_dimob/hmmer-3.1b1/bin:" + self.config.SOFTWARE_DIR + "/bioinfo/metaGenomic/pplacer-Linux-v1.1.alpha19"
        self.set_environ(PATH=self.path)
        self.checkm = "/program/Python/bin/"
        self.perl = "/program/perl/perls/perl-5.24.0/bin/perl"
        self.perl_script = self.config.PACKAGE_DIR + "/metagbin/checkm_stat.pl"
        self.bin_dir =self.option('bin_dir').prop['path']
        self.summary =self.work_dir + '/checkm.summary.txt'
        self.result = self.work_dir + '/result'
        self.taxon = self.config.SOFTWARE_DIR + "/bioinfo/metaGenomic/CheckM-master/checkm.taxon.xls"
        self.bin_name = os.listdir(self.bin_dir)[0].split('.')[0]

    def run_checkm(self):
        if os.path.exists(self.result):
            shutil.rmtree(self.result)
        cmd = "{}checkm lineage_wf -x fa -r -t 16 --pplacer_threads 16 -f {} --tab_table --tmpdir {} {} {}".format(self.checkm,self.summary, "./",self.bin_dir,self.result)
        self.logger.info(cmd)
        self.logger.info("开始运行run_checkm")
        command = self.add_command("run_checkm", cmd)
        command.run()
        self.wait(command)
        if command.return_code == 0:
            self.logger.info("运行run_checkm完成")
        else:
            self.set_error("运行run_checkm运行出错!")

    def run_checkm_stat(self):
        cmd = "{} {} {} {} {} {}".format(self.perl,self.perl_script,self.summary,self.result + '/storage/bin_stats_ext.tsv',self.taxon,self.bin_name)
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
        for i in [self.bin_name + '.marker.xls',self.bin_name + '.bin.summary.xls']:
            if os.path.exists(self.output_dir + "/" + i):
                os.remove(self.output_dir + "/" + i)
            os.link(self.work_dir + "/" + i,self.output_dir + "/" + i)

    def run(self):
        """
        运行
        """
        super(CheckmBinTool, self).run()
        self.run_checkm()
        self.run_checkm_stat()
        self.set_output()
        self.end()