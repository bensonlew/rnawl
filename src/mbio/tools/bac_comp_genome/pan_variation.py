# -*- coding: utf-8 -*-
# __author__ = 'qingchen.zhang'
# version 1.0
# last_modify: 20119.10.10

import os
import re, shutil
import pandas as pd
from biocluster.agent import Agent
from biocluster.tool import Tool
from biocluster.core.exceptions import OptionError
from mbio.packages.bac_comp_genome.common_function import link_dir


class PanVariationAgent(Agent):
    """
    细菌比较基因组变异分析，原PGAP代码没有改动
    如果是pgap预测的结果，就不需要再重新整理了，直接用pgap软件上一步聚类结果的输入文件就可以了
    然后用整理后的结果
    """
    def __init__(self, parent):
        super(PanVariationAgent, self).__init__(parent)
        options = [
            {"name": "infile_dir", "type": "infile", "format": "bac_comp_genome.input_dir"},  # 输入原始文件夹.faa和.fna文件
            {"name": "cluster", "type": "infile", "format": "sequence.profile_table"},  # 输入聚类的cluster结果整理的二维结果表
        ]
        self.add_option(options)

    def check_options(self):
        """
        重写参数检查
        :return:
        """
        if not self.option("infile_dir").is_set:
            raise OptionError("必须设置参数infile_dir文件夹!")
        if not self.option("cluster").is_set:
            raise OptionError("必须设置聚类结果二维表")

    def set_resource(self):
        self._cpu = 3
        self._memory = '30G'

    def end(self):
        super(PanVariationAgent, self).end()

class PanVariationTool(Tool):
    def __init__(self, config):
        super(PanVariationTool, self).__init__(config)
        self.path = self.config.SOFTWARE_DIR + "/bioinfo/compare_genome/software/PGAP-1.2.1:" + self.config.SOFTWARE_DIR + "/bioinfo/compare_genome/software/mcl/bin:" + self.config.SOFTWARE_DIR + \
                    "/bioinfo/align/ncbi-blast-2.3.0+/bin:" + self.config.SOFTWARE_DIR + "bioinfo/phylogenetic/standard-RAxML-master:" + self.config.SOFTWARE_DIR + "/bioinfo/compare_genome/software/phylip-3.697/exe:"
        self.perl5path = self.config.SOFTWARE_DIR + "/bioinfo/compare_genome/software/PGAP-1.2.1"
        self.set_environ(PATH=self.path, PERL5LIB=self.perl5path)
        #self.strains = self.option("strains")
        self.fasta = self.option("infile_dir").prop['path']
        self.perl_path = "/program/perl/perls/perl-5.24.0/bin/perl"
        self.perl_script = self.config.SOFTWARE_DIR + "/bioinfo/compare_genome/software/PGAP-1.2.1/PGAP.pl"
        self.out = self.work_dir + '/out_result'
        self.python = "/program/Python/bin/python"
        self.format_convert = self.config.PACKAGE_DIR + "/bac_comp_genome/format_convert.py"

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
        cmd = "{} {} -i {} -out {} -o {} -t 3".format(self.python, self.format_convert, self.fasta, self.fasta_dir, self.input_dir)
        self.logger.info(cmd)
        command = self.add_command("format_convert", cmd).run()
        self.wait(command)
        if command.return_code == 0:
            self.logger.info("格式转换运行完成！")
        else:
            self.set_error("格式转换运行出错!")

    def run_variation(self):
        """
        用pgap软件进行变异分析
        :return:
        """
        self.logger.info("开始用PGAP软件进行变异分析")
        self.get_strains()
        input_cluster = self.option("cluster").prop['path']
        cluster_name = '1.Orthologs_Cluster.txt'
        cluster_path = os.path.join(self.out, cluster_name)
        if os.path.exists(self.out):
            shutil.rmtree(self.out)
        os.mkdir(self.out)
        if os.path.exists(cluster_path):
            os.remove(cluster_path)
        os.link(input_cluster, cluster_path)
        new_cluster_path = os.path.join(self.out, 'new_1.Orthologs_Cluster.txt')
        self.remove_columns(cluster_path, new_cluster_path)
        os.remove(cluster_path)
        os.system("mv %s %s"%(new_cluster_path, cluster_path))
        cmd = "{} {} --strains {} --input {} --output {} --variation --thread 3".format(self.perl_path, self.perl_script, self.strains, self.input_dir, self.out)
        self.logger.info(cmd)
        command = self.add_command("run_variation", cmd).run()
        self.wait(command)
        if command.return_code == 0:
            self.fix_rownames()
            self.logger.info("run_variation运行完成！")
        else:
            self.set_error("run_variation运行完成运行出错!")

    def fix_rownames(self):
        """
        对变异分析运行的结果进行改名
        规则是将第一列的数字前加上CLUSTER
        :return:
        """
        variation_path = os.path.join(self.out, '3.CDS.variation.txt')
        new_variation_path = os.path.join(self.out, 'new_3.CDS.variation.txt')
        analysis_path = os.path.join(self.out, '3.CDS.variation.analysis.txt')
        new_analysis_path = os.path.join(self.out, 'new_3.CDS.variation.analysis.txt')
        with open(variation_path, "r") as fr, open(new_variation_path, 'w') as w:
            lines = fr.readlines()
            w.write(lines[0])
            for line in lines[1:]:
                line = line.strip().split("\t")
                cluster_name = "CLUSTER"+str(line[0])
                line[0] = cluster_name
                w.write("\t".join(line) +"\n")
        with open(analysis_path, "r") as fa, open(new_analysis_path, 'w') as wx:
            n_lines = fa.readlines()
            wx.write(n_lines[0])
            for lin in n_lines[1:]:
                lin = lin.strip().split("\t")
                cluster_name = "CLUSTER"+str(lin[0])
                lin[0] = cluster_name
                wx.write("\t".join(lin) +"\n")

    def get_strains(self):
        """
        根据传入的文件夹内的文件生成对应的样本名称
        :return:
        """
        cluster_path = self.option("cluster").prop["path"]
        sample_list = []
        if os.path.exists(cluster_path):
            with open(cluster_path, "r") as f:
                line = f.readline()
                sample_list = line.strip().split("\t")[3:]
        self.strains = "+".join(sample_list)

    def remove_columns(self, cluster, new_cluster):
        """
        主要是删除第2、3列的统计信息，并将第一列的名称改写，便于进行变异分析
        修改bug统计错误导致的原因
        :return:
        """
        with open(cluster, 'r') as f, open(new_cluster, 'w') as w:
            lines = f.readlines()
            header_list = lines[0].strip().split("\t")
            header_firt = []
            header_firt.append(header_list[0])
            new_header_list = header_firt + header_list[3:]
            new_header_name = "\t".join(new_header_list)
            w.write(new_header_name + "\n")
            for line in lines[1:]:
                line = line.strip().split("\t")
                name = line[0].strip("CLUSTER")
                cluster_name = []
                cluster_name.append(str(name))
                sample_list = line[3:]
                total_list = cluster_name + sample_list
                total_sample = "\t".join(total_list)
                w.write(total_sample+ "\n")
        # data = pd.read_table(cluster, sep='\t',header=0)
        # list = []
        # with open(cluster, 'r') as f:
        #     line = f.readline()
        #     line = line.strip().split("\t")
        #     list = line
        # list.pop(1)
        # list.pop(1)
        # new_data = data[list]
        # new_data.to_csv(new_cluster, sep="\t",index=0)

    def set_output(self):
        """
        设置结果文件目录,对结果文件名称改名字
        :return:
        """
        variation_path = os.path.join(self.output_dir, 'CDS_variation.xls')
        analysis_path = os.path.join(self.output_dir, 'CDS_variation_analysis.xls')
        if os.path.exists(variation_path):
            os.remove(variation_path)
        if os.path.exists(analysis_path):
            os.remove(analysis_path)
        os.link(os.path.join(self.out, 'new_3.CDS.variation.txt'), variation_path)
        os.link(os.path.join(self.out, 'new_3.CDS.variation.analysis.txt'), analysis_path)

    def run(self):
        """
        运行
        :return:
        """
        self.logger.info("开始运行tool")
        super(PanVariationTool, self).run()
        self.format_covert()
        self.run_variation()
        self.set_output()
        self.end()
