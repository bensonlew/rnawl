# -*- coding: utf-8 -*-
# __author__ = 'qingchen.zhang'
# version 1.0
# last_modify: 20119.10.10

import os
import random
import re, shutil
import pandas as pd
from biocluster.agent import Agent
from biocluster.tool import Tool
from biocluster.core.exceptions import OptionError
from mbio.packages.bac_comp_genome.common_function import link_dir

class PanGenomeAgent(Agent):
    """
    细菌比较基因组计算泛基因组大小公式，修改改写PGAP这部分的代码运算方式，主要是加入取样方法，使得分析更加的合理，而且要加入新基因的计算方法
    """
    def __init__(self, parent):
        super(PanGenomeAgent, self).__init__(parent)
        options = [
            {"name": "infile_dir", "type": "infile", "format": "bac_comp_genome.input_dir"},  # 输入原始文件夹.faa和.fna文件
            {"name": "cluster", "type": "infile", "format": "sequence.profile_table"},  # 输入聚类的cluster结果整理的二维结果表
            {"name": "pangenome", "type": "string", "default": "homologus"}, ## 计算泛基因组大小的公式的方法
            {"name": "cal_method", "type": "string", "default": "core_Tettelin"}, ###计算方法
            {"name": "cal_newgene", "type": "string", "default": "average"}, ##计算newgene的方法
        ]
        self.add_option(options)

    def check_options(self):
        """
        重写参数检查
        :return:
        """
        # if not self.option("infile_dir").is_set:
        #     raise OptionError("必须设置参数infile_dir文件夹!")
        if not self.option("cluster").is_set:
            raise OptionError("必须设置聚类结果二维表")

    def set_resource(self):
        self._cpu = 2
        self._memory = '20G'

    def end(self):
        super(PanGenomeAgent, self).end()

class PanGenomeTool(Tool):
    def __init__(self, config):
        super(PanGenomeTool, self).__init__(config)
        self.path = self.config.SOFTWARE_DIR + "/bioinfo/compare_genome/software/PGAP-1.2.1:" + self.config.SOFTWARE_DIR + "/bioinfo/compare_genome/software/mcl/bin:" + self.config.SOFTWARE_DIR + \
                    "/bioinfo/align/ncbi-blast-2.3.0+/bin:" + self.config.SOFTWARE_DIR + "bioinfo/phylogenetic/standard-RAxML-master:" + self.config.SOFTWARE_DIR + "/bioinfo/compare_genome/software/phylip-3.697/exe:"
        self.perl5path = self.config.SOFTWARE_DIR + "/bioinfo/compare_genome/software/PGAP-1.2.1"
        self.set_environ(PATH=self.path, PERL5LIB=self.perl5path)
        # if self.option("pangenome") == 'pgap':
        self.fasta = self.option("infile_dir").prop['path']
        self.perl_path = "/miniconda2/bin/perl"
        self.perl_script = self.config.SOFTWARE_DIR + "/bioinfo/compare_genome/software/PGAP-1.2.1/PGAP.pl"
        self.out = self.work_dir + '/out_result'
        self.python = "/miniconda2/bin/python"
        self.format_convert = self.config.PACKAGE_DIR + "/bac_comp_genome/format_convert.py"
        self.calculate_pangenome = self.config.PACKAGE_DIR + "/bac_comp_genome/calculate_pangenome.py"
        self.newgene = self.config.PACKAGE_DIR + "/bac_comp_genome/calculate_newgene.py"
        self.path = self.config.SOFTWARE_DIR + "/bioinfo/compare_genome/software/get_homologues-master:" + self.config.SOFTWARE_DIR + "/program/R-3.3.1/bin:"
        self.perl5path = self.config.SOFTWARE_DIR + "/miniconda2/bin/perl"
        self.library = self.config.SOFTWARE_DIR + "/program/R-3.3.1/lib64/R/lib"
        self.r_home = self.config.SOFTWARE_DIR + "/program/R-3.3.1/lib64/R"
        self.c_include = self.config.SOFTWARE_DIR + "/program/R-3.3.1/lib64/R/include"
        self.set_environ(PATH=self.path, PERL5LIB=self.perl5path)
        self.core_plot = self.config.SOFTWARE_DIR + "/bioinfo/compare_genome/software/get_homologues-master/plot_pancore_matrix.pl"

    def format_covert(self):
        """
        对输入文件进行标准化处理，为序列名称加样本名称并以“|”分割，为变异分析整理成pgap的输入格式.pep .nuc .function
        .function文件为功能分析的文件，现转了一个只有标题的文件
        :return:
        """

        self.logger.info("开始进行格式转换")
        self.fasta_dir = self.work_dir + '/fasta_dir'
        if os.path.exists(self.fasta_dir):
            shutil.rmtree(self.fasta_dir)
        os.mkdir(self.fasta_dir)
        self.input_dir = self.work_dir + '/input_dir'
        if os.path.exists(self.input_dir):
            shutil.rmtree(self.input_dir)
        os.mkdir(self.input_dir)
        cmd = "{} {} -i {} -out {} -o {} -t 2".format(self.python, self.format_convert, self.fasta, self.fasta_dir, self.input_dir)
        self.logger.info(cmd)
        command = self.add_command("format_convert", cmd).run()
        self.wait(command)
        if command.return_code == 0:
            self.logger.info("格式转换运行完成！")
        else:
            self.set_error("格式转换运行出错!")

    def run_pangenome(self):
        """
        运行PGAP进行泛基因组大小的计算
        :return:
        """
        sample_list = self.get_strains()
        sample_num = len(sample_list)
        if sample_num >= 60:
            sampleSize = 50
        else:
            sampleSize = 100
        input_cluster = self.option("cluster").prop['path']
        cluster_name = '1.Orthologs_Cluster.txt'
        cluster_path = os.path.join(self.out, cluster_name)
        if os.path.exists(cluster_path):
            os.remove(cluster_path)
        if not os.path.exists(self.out):
            os.mkdir(self.out)
        os.link(input_cluster, cluster_path)
        new_cluster_path = os.path.join(self.out, 'new_1.Orthologs_Cluster.txt')
        self.remove_columns(cluster_path, new_cluster_path)
        os.remove(cluster_path)
        os.system("mv %s %s"%(new_cluster_path, cluster_path))
        cmd = "{} {} --strains {} --input {} --output {} --pangenome --sampleSize {} --thread 2 ".format(self.perl_path, self.perl_script, self.strains, self.input_dir, self.out, sampleSize)

        command = self.add_command("run_pangenome", cmd).run()
        self.wait(command)
        if command.return_code == 0:
            self.logger.info("run_pangenome运行完成！")
        else:
            self.set_error("run_pangenome运行完成运行出错!")

    def get_strains(self):
        """
        根据传入的文件夹内的文件生成对应的样本名称
        :return:
        """
        all_files = os.listdir(self.fasta_dir)
        sample_list = []
        for file in all_files:
            if re.search(r".pep", file):
                sample = file.strip("_CDS.pep")
                sample_list.append(sample)
        self.strains = "+".join(sample_list)
        return sample_list

    def remove_columns(self, cluster, new_cluster):
        """
        主要是删除第2、3列的统计信息，便于进行变异分析
        :return:
        """
        data = pd.read_table(cluster, sep='\t',header=0)
        list = []
        with open(cluster, 'r') as f:
            line = f.readline()
            line = line.strip().split("\t")
            list = line
        list.pop(1)
        list.pop(1)
        new_data = data[list]
        new_data.to_csv(new_cluster, sep="\t",index=0)

    def calculate_genome(self):
        """
        主要根据输入的cluster表进行整理和计算出get_homologus的公式
        :return:
        """
        if os.path.exists(self.out):
            shutil.rmtree(self.out)
        os.mkdir(self.out)
        cmd = "{} {} -i {} -out {}".format(self.python, self.calculate_pangenome, self.option("cluster").prop['path'], self.out)
        self.logger.info(cmd)
        command = self.add_command("calculate_pangenome", cmd).run()
        self.wait(command)
        if command.return_code == 0:
            self.logger.info("计算core和pan运行完成！")
        else:
            self.set_error("计算pan和core运行出错!")

    def calculate_newgene(self):
        """
        计算newgene
        :return:
        """
        outfile = os.path.join(self.out, "new_gene_cluster.xls.log")
        if os.path.exists(outfile):
            os.remove(outfile)
        infile_path = os.path.join(self.out, "new_gene_cluster.xls")
        cmd = "{} {} -i {} -o {} -m {}".format(self.python, self.newgene, infile_path, outfile, self.option('cal_newgene'))
        self.logger.info(cmd)
        command = self.add_command("calculate_newgene", cmd).run()
        self.wait(command)
        if command.return_code == 0:
            self.logger.info("计算newgene运行完成！")
        else:
            self.set_error("计算newgene运行出错!")

    def plot_core(self, method):
        """
        core和pan计算公式
        :return:
        """
        if method in ['pan']:
            infile_path = os.path.join(self.out, 'pan_cluster.xls')
        elif method in ["core_Tettelin", "core_Willenbrock"]:
            infile_path = os.path.join(self.out, 'core_cluster.xls')
        else:
            infile_path = os.path.join(self.out, 'new_gene_cluster.xls')
        cmd = "{} {} -i {} -f {}".format(self.perl_path, self.core_plot, infile_path, method)
        self.logger.info(cmd)
        command_name = "plot_core" + method.lower()
        command = self.add_command(command_name, cmd).run()
        self.wait(command)
        if command.return_code == 0:
            self.logger.info("计算公式运行完成！")
        else:
            self.set_error("计算公式运行出错!")

    def set_output(self):
        """
        设置结果文件目录
        :return:
        """
        self.logger.info("正在设置结果文件目录")
        for file in os.listdir(self.out):
            file_path = os.path.join(self.out, file)
            newfile = os.path.join(self.output_dir, file)
            if os.path.exists(newfile):
                os.remove(newfile)
            if file.endswith(".log"):
                os.link(file_path, os.path.join(self.output_dir, file))
            elif file.endswith(".xls"):
                os.link(file_path, os.path.join(self.output_dir, file))
        genome_path = os.path.join(self.output_dir, "pgap_formula.xls")
        if os.path.exists(genome_path):
            os.remove(genome_path)
        if os.path.exists(os.path.join(self.out, "2.PanGenome.Profile.txt")):
            os.link(os.path.join(self.out, "2.PanGenome.Profile.txt"), genome_path)
        self.logger.info("设置结果目录文件完成")

    def run(self):
        super(PanGenomeTool, self).run()
        if self.option("pangenome") == 'pgap':
            self.format_covert()
            self.run_pangenome()
            # self.calculate_genome()
            # self.calculate_newgene()
        else:
            self.calculate_genome()
            for i in [self.option("cal_method"), 'pan']:
                self.plot_core(i)
            self.calculate_newgene()
            # self.format_covert()
            # self.run_pangenome()
        self.set_output()
        self.end()



