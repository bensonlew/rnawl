#!/usr/bin/env python
# -*- coding: utf-8 -*-
from biocluster.agent import Agent
from biocluster.tool import Tool
from biocluster.core.exceptions import OptionError
import os
import shutil
from mbio.packages.metagbin.common_function import link_dir


class ResfinderAgent(Agent):
    """
    用于resfinder进行耐药基因的预测
    author: qingchen.zhang
    """

    def __init__(self, parent):
        super(ResfinderAgent, self).__init__(parent)
        options = [
            {"name": "gene_gff", "type": "infile", "format": "gene_structure.gff3"},  #样品gff文件
            {"name": "gene_fa", "type": "infile", "format": "sequence.fasta"},  #样品gene序列
            {"name": "min_cov", "type": "float", "default": 0.6},# 最小的coverage
            {"name": "min_iden", "type": "float", "default": 0.8},## 最小的identity
            {"name": "sample_name", "type": "string"}, ## 样品名称
            {"name": "resfinder_database", "type": "string", "default": "true"},## 有三个数据库，默认做resfinder
            {"name": "pointfinder_database", "type": "string", "default": "false"},
            {"name": "disinfinder_database", "type": "string", "default": "true"},
            {"name": "species_name", "type": "string", "default": ""}, ## 点突变物种名称
        ]
        self.add_option(options)
        self.list =[]


    def check_options(self):
        """
        检测参数是否正确
        """
        if not self.option('gene_fa').is_set:
            raise OptionError("请设置基因组基因序列不存在！")


    def set_resource(self):
        """
        所需资源
        """
        self._cpu = 1
        self._memory = '20G'

    def end(self):
        super(ResfinderAgent, self).end()


class ResfinderTool(Tool):
    """
    version 1.0
    """
    def __init__(self, config):
        super(ResfinderTool, self).__init__(config)
        self.soft_path = self.config.SOFTWARE_DIR + "/program/Python35/bin:" + self.config.SOFTWARE_DIR + \
                         "/gcc/5.1.0/bin:"
        self.soft_lib = self.config.SOFTWARE_DIR + "/program/Python35/lib:" + self.config.SOFTWARE_DIR + "/bioinfo/Genomic/Sofware/resfinder:" + self.config.SOFTWARE_DIR + "/gcc/5.1.0/lib64:" + self.config.SOFTWARE_DIR + "/program/Python/lib:"
                        # self.config.SOFTWARE_DIR + "/bioinfo/Genomic/Sofware/resfinder/.git/"
        self.set_environ(PATH=self.soft_path, LD_LIBRARY_PATH=self.soft_lib)
        # self.set_environ(GIT_PYTHON_REFRESH="quiet")
        self.seq = self.option('gene_fa').prop['path']
        self.gff = self.option('gene_gff').prop['path']
        self.blastn = self.config.SOFTWARE_DIR + "/bioinfo/align/ncbi-blast-2.3.0+/bin/blastn"
        self.resfinder_database = self.config.SOFTWARE_DIR + "/bioinfo/Genomic/Sofware/resfinder/db_resfinder"
        self.pointfinder_database = self.config.SOFTWARE_DIR + "/bioinfo/Genomic/Sofware/resfinder/db_pointfinder"
        self.disinfinder_database = self.config.SOFTWARE_DIR + "/bioinfo/Genomic/Sofware/resfinder/db_disinfinder"
        self.python3 = "/program/Python35/bin/python3"
        self.scripts = self.config.SOFTWARE_DIR + "/bioinfo/Genomic/Sofware/resfinder/run_resfinder.py"
        self.resfinder_pheno = self.config.SOFTWARE_DIR + "/bioinfo/Genomic/Sofware/resfinder/db_resfinder/phenotypes.txt"
        self.disinfinder_pheno = self.config.SOFTWARE_DIR + "/bioinfo/Genomic/Sofware/resfinder/db_disinfinder/phenotypes.txt"
        self.python = "/miniconda2/bin/python"
        self.stat = self.config.PACKAGE_DIR + "/bacgenome/resfinder_stat.py"


    def run_resfinder(self):
        """
        运行resfinder软件进行预测和比对
        """
        output = self.work_dir + "/out_result"
        if os.path.exists(output):
            shutil.rmtree(output)
        os.mkdir(output)
        cmd = "{} {} -ifa {} -o {} --acquired --blastPath {}".format(self.python3, self.scripts, self.seq, output, self.blastn)
        if self.option("resfinder_database") in ['true']:
            cmd += " --min_cov {} -t {}".format(self.option("min_cov"), self.option("min_iden"))
            cmd += " --db_path_res {}".format(self.resfinder_database)
        if self.option("pointfinder_database") in ['true']:
            cmd += " --point "
            cmd += " --min_cov_point {} --threshold_point {}".format(self.option("min_cov"), self.option("min_iden"))
            cmd += " --db_path_point {}".format(self.pointfinder_database)
        if self.option("disinfinder_database") in ['true']:
            cmd += " --disinfectant "
            cmd += " --db_path_disinf {}".format(self.disinfinder_database)
        if self.option("species_name") not in ['']:
            cmd += " -s {} ".format(self.option("species_name"))

        self.logger.info(cmd)
        self.logger.info("开始运行resfinder")
        command = self.add_command("run_resfinder", cmd)
        command.run()
        self.wait(command)
        if command.return_code == 0:
            self.logger.info("运行run_resfinder完成")
        else:
            self.set_error("运行run_resfinder运行出错!")

    def run_stat(self):
        """
        对resfinder的结果进行整理和统计
        """
        resfinder_results = self.work_dir + "/out_result/ResFinder_results_tab.txt"
        stat = self.work_dir + "/stat"
        if os.path.exists(stat):
            shutil.rmtree(stat)
        cmd = "{} {} -i {} -s {} -o {} -n {} -gff {} -t {}".format(self.python, self.stat, resfinder_results,self.resfinder_pheno, stat, self.option("sample_name"), self.gff, "resfinder")
        self.logger.info(cmd)
        self.logger.info("开始运行resfinder")
        command = self.add_command("run_stat", cmd)
        command.run()
        self.wait(command)
        if command.return_code == 0:
            self.logger.info("运行run_stat完成")
        else:
            self.set_error("运行run_stat运行出错!")

    def run_disinfinder(self):
        """
        对disinfinder结果进行整理
        :return:
        """
        disinfinder_results = self.work_dir + "/out_result/DisinFinder_results_tab.txt"
        stat = self.work_dir + "/out_stat"
        if os.path.exists(stat):
            shutil.rmtree(stat)
        cmd = "{} {} -i {} -s {} -o {} -n {} -gff {} -t {}".format(self.python, self.stat, disinfinder_results,
                                                             self.disinfinder_pheno, stat, self.option("sample_name"),
                                                             self.gff, "disinfinder")
        self.logger.info(cmd)
        self.logger.info("开始运行resfinder")
        command = self.add_command("run_stat2", cmd)
        command.run()
        self.wait(command)
        if command.return_code == 0:
            self.logger.info("运行run_stat完成")
        else:
            self.set_error("运行run_stat运行出错!")

    def merge_stat(self):
        """
        合并统计文件
        """
        if os.path.exists(self.work_dir + "/out_stat") and os.path.exists(self.work_dir + "/stat"):
            a1 = self.work_dir + "/out_stat"
            a2 = self.work_dir + "/stat"
            stat1 = os.path.join(a1, self.option("sample_name") + ".stat.xls")
            stat2 = os.path.join(a2, self.option("sample_name") + ".stat.xls")
            with open(self.work_dir + "/" + self.option("sample_name") + ".stat.xls", 'w') as w, open(stat1, 'r') as f1,open(stat2, 'r') as f2:
                w.write("Sample Name\tGene No.\n")
                f1.readline()
                f2.readline()
                f1_line = f1.readline().strip().split("\t")[1]
                f2_line = f2.readline().strip().split("\t")[1]
                all_gene_number = int(f1_line) + int(f2_line)
                w.write("{}\t{}\n".format(self.option("sample_name"), all_gene_number))


    def set_output(self):
        """
        设置结果文件目录
        """
        self.logger.info("开始生成结果文件目录")
        path =self.output_dir + '/' + self.option("sample_name")
        if os.path.exists(path):
            shutil.rmtree(path)
        if os.path.exists(self.work_dir + "/stat"):
            link_dir(self.work_dir + "/stat", path)
        if os.path.exists(self.work_dir + "/out_stat"):
            link_dir(self.work_dir + "/out_stat", path)
        if os.path.exists(self.work_dir + "/" + self.option("sample_name") + ".stat.xls"):
            if os.path.exists(path + "/" + self.option("sample_name") + ".stat.xls"):
                os.remove(path + "/" + self.option("sample_name") + ".stat.xls")
            os.link(self.work_dir + "/" + self.option("sample_name") + ".stat.xls", path + "/" + self.option("sample_name") + ".stat.xls")


    def run(self):
        """
        运行
        """
        super(ResfinderTool, self).run()
        self.run_resfinder()
        self.run_stat()
        self.run_disinfinder()
        self.merge_stat()
        self.set_output()
        self.end()
