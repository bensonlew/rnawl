#!/usr/bin/env python
# -*- coding: utf-8 -*-
from biocluster.agent import Agent
from biocluster.tool import Tool
from biocluster.core.exceptions import OptionError
import os,re,shutil
from mbio.packages.bac_comp_genome.common_function import link_file,link_dir
import subprocess


class BacGbkAgent(Agent):
    """
    用于细菌的gbk文件生成，再进行次级代谢产物分析,基因组岛分析
    version 1.0
    author: gaohao
    last_modify: 2019.09.28
    """

    def __init__(self, parent):
        super(BacGbkAgent, self).__init__(parent)
        options = [
            {"name": "gene_gff", "type": "infile", "format": "gene_structure.gff3"},  #
            {"name": "rrna_gff", "type": "infile", "format": "gene_structure.gff3"},  #
            {"name": "trna_gff", "type": "infile", "format": "gene_structure.gff3"},  #
            {"name": "pro_fa", "type": "infile", "format": "sequence.fasta"},  #
            {"name": "genome_fa", "type": "infile", "format": "sequence.fasta"},  #
            {"name": "sample_name", "type": "string"},
            {"name": "path", "type": "string"},
            {"name": "type", "type": "string", "default": "genome"},
            {"name": "gbk_dir", "type": "outfile", "format": "gene_structure.gbk_dir"},  # 输出文件
            {"name": "all_gbk", "type": "outfile", "format": "gene_structure.gbk"},  # 输出文件
        ]
        self.add_option(options)

    def check_options(self):
        """
        检测参数是否正确
        """
        if not self.option('gene_gff').is_set:
            raise OptionError("请设置基因组基因预测基因gff文件！")
        if not self.option('rrna_gff').is_set:
            raise OptionError("请设置基因组rRNA预测基因gff文件！")
        if not self.option('trna_gff').is_set:
            raise OptionError("请设置基因组tRNA预测基因gff文件！")
        if not self.option('pro_fa').is_set:
            raise OptionError("请设置基因组基因蛋白序列文件！")
        if not self.option('genome_fa').is_set:
            raise OptionError("请设置基因组组装序列文件！")


    def set_resource(self):
        """
        所需资源
        """
        self._cpu = 3
        self._memory = '10G'

    def end(self):
        super(BacGbkAgent, self).end()


class BacGbkTool(Tool):
    """
    version 1.0
    """
    def __init__(self, config):
        super(BacGbkTool, self).__init__(config)
        self.genegff =self.option('gene_gff').prop['path']
        self.rrnagff = self.option('rrna_gff').prop['path']
        self.trnagff = self.option('trna_gff').prop['path']
        self.genomefa = self.option('genome_fa').prop['path']
        self.pro = self.option('pro_fa').prop['path']
        self.sample = self.option("sample_name")
        self.perl_path = "/program/perl-5.24.0/bin/perl"
        self.perl_script = self.config.PACKAGE_DIR + "/bac_comp_genome/"

    def run_bio_gbk(self):
        if os.path.exists(self.work_dir + "/gbk"):
            shutil.rmtree(self.work_dir + "/gbk")
        cmd = "{} {}GBK_generation.pl {} {} {} {} {}".format(self.perl_path, self.perl_script,
                                                                            self.genegff, self.trnagff, self.rrnagff,
                                                                            self.pro, self.genomefa)
        self.logger.info(cmd)
        self.logger.info("开始运行run_bio_gbk")
        command = self.add_command("run_bio_gbk", cmd)
        command.run()
        self.wait(command)
        if command.return_code == 0:
            self.logger.info("运行run_bio_gbk完成")
        else:
            self.set_error("运行run_bio_gbk运行出错!")

    def run_bio_gbk2(self):
        if os.path.exists(self.work_dir + "/gbk"):
            shutil.rmtree(self.work_dir + "/gbk")
        cmd = "{} {}GBK_generation2.pl {} {} {} {} {} {}".format(self.perl_path, self.perl_script,
                                                             self.genegff, self.trnagff, self.rrnagff,
                                                             self.pro, self.genomefa, self.option("path"))
        self.logger.info(cmd)
        self.logger.info("开始运行run_genome_gbk")
        command = self.add_command("run_genome_gbk", cmd)
        command.run()
        self.wait(command)
        if command.return_code == 0:
            self.logger.info("运行run_genome_gbk完成")
        else:
            self.set_error("运行run_genome_gbk运行出错!")

    def cat_gbk(self):
        cmd = "cat {} >{}".format(self.work_dir + "/gbk/*", self.work_dir + "/" + self.sample + '.gbk')
        try:
            self.logger.info(cmd)
            subprocess.check_output(cmd, shell=True)
            self.logger.info("cat_gbk运行完成")
        except subprocess.CalledProcessError:
            self.set_error("cat_gbk运行出错")

    def set_output(self):
        self.logger.info("set output")
        link_file(self.work_dir + "/" + self.sample + '.gbk', self.output_dir + "/" + self.sample + '.gbk')
        self.option('all_gbk').set_path(self.output_dir + "/" + self.sample + '.gbk')
        link_dir(self.work_dir + "/gbk", self.output_dir + '/' + self.sample)
        self.option('gbk_dir').set_path(self.output_dir + '/' + self.sample)

    def run(self):
        """
        运行
        """
        super(BacGbkTool, self).run()
        if self.option('type') in ['genome', 'Genome']:
            self.run_bio_gbk2()
        else:
            self.run_bio_gbk()
        self.logger.info("aaaa")
        self.cat_gbk()
        self.logger.info("bbbb")
        self.set_output()
        self.end()