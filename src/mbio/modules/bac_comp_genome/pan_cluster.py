# -*- coding: utf-8 -*-
#__author__ = 'qingchen.zhang'
import os
import re
import time
from biocluster.core.exceptions import OptionError
from biocluster.module import Module
from collections import defaultdict
import shutil
import gevent
from mbio.packages.bac_comp_genome.common_function import link_dir


class PanClusterModule(Module):
    def __init__(self, work_id):
        """
        细菌比较基因组泛基因组分析--聚类分析用的module
        然后进行pangenome大小的计算和变异分析
        :return:
        """
        super(PanClusterModule, self).__init__(work_id)
        options = [
            {"name": "fasta_dir", "type": "infile", "format": "bac_comp_genome.input_dir"},  # 输入核酸fasta的dir文件,里面包含faa文件和fna文件
            {"name": "identity", "type": "float", "default": 0.5},  ##给出cdhit的参数identity
            {"name": "coverage", "type": "float", "default": 0.0},  # 给出cdhit的参数coverage
            {"name": "method", "type": "string", "default":'roary'}, #输入方法类型{'orthofinder', 'orthomcl', 'get_homologus', 'pgap', 'roary'}
            {"name": "cluster_method", "type": "string", "default": "cdhit+blast+mcl"}, #cdhit usearch mmseq blast+mcl
            {"name": "inflation", "type": "float", "default": 1.5},  # Markov Inflation Index
            {"name": "out", "type": "outfile", "format": "sequence.profile_table"},  # 输出二维表格
            {"name": "pgap_mehod", "type": "string", "default": "GF"},  # GF for GeneFamily method,  and MP for MultiParanoid method
            {"name": "homologus_mehod", "type": "string", "default": "OMCL"},  # 三种方法BDBH、OMCL、OCOG
            # {"name": "anno_file", "type": "infile", "format": "sequence.profile_table"},  # 物种注释总览表
            {"name": "score", "type": "int", "default": 40},#Minimum score in blast
            {"name": "evalue", "type": "string", "default": "1e-10"},#Maximal E-value in blastall
            {"name": "local", "type": "float", "default": 0.25},#Minimum local alignment overlap in MP method
            {"name": "global", "type": "float", "default": 0.5},#Minimum global alignment overlap in MP method
        ]
        self.add_option(options)
        self.convert_format = self.add_tool("bac_comp_genome.add_sample")
        self.cluster_tools = []

    def check_options(self):
        """
        进行参数二次检查
        :return:
        """
        if not self.option("fasta_dir").is_set:
            raise OptionError("必须设置参数fasta_dir")
        if not self.option("cluster_method"):
            raise OptionError("必须设置聚类方法cluster_method")
        return True

    def convert_data(self):
        """
        将序列加上样本名称，如果为cdhit、mmseq和usearch则需要合并
        :return:
        """
        self.logger.info("正在用进行格式转换")
        opts=({
            "fasta_dir": self.option("fasta_dir"),
            })
        if self.option("cluster_method") in ["cdhit", "usearch", "mmseq"]:
            opts["is_merge"] = 1
        else:
            opts["is_merge"] = 0
        self.convert_format.set_options(opts)
        self.convert_format.run()

    def run_mmseq(self):
        """
        用mmseq软件进行聚类
        :return:
        """
        self.logger.info("正在用mmseq进行聚类")
        self.mmseq = self.add_tool("metagbin.mmseq_seq")
        query_fasta = self.convert_format.option('out')
        self.mmseq.set_options({
            "contig_fa": query_fasta,
            "cdhit_identity": self.option("identity"),
            "cdhit_coverage": self.option("coverage"),
            })
        self.mmseq.on('end', self.set_output)
        self.mmseq.run()
        self.cluster_tools.append(self.mmseq)

    def run_cdhit(self):
        """
        用cdhit软件进行聚类
        :return:
        """
        self.logger.info("正在用cdhit进行聚类")
        self.compare_single = self.add_tool("bac_comp_genome.pan_cdhit")
        query_fasta = self.convert_format.option('out')
        self.compare_single.set_options({
            "query": query_fasta,
            "identity": self.option("identity"),
            "coverage": self.option("coverage"),
            "ana_type": "prot"
            })
        self.compare_single.on('end', self.set_output)
        self.compare_single.run()
        self.cluster_tools.append(self.compare_single)

    def run_usearch(self):
        """
        用usearch进行聚类
        :return:
        """
        self.logger.info("正在用usearch进行聚类")
        self.usearch = self.add_tool("bac_comp_genome.pan_usearch")
        query_fasta = self.convert_format.option('out')
        self.usearch.set_options({
            "fasta": query_fasta,
            "identity": self.option("identity"),
            })
        self.usearch.on('end', self.set_output)
        self.usearch.run()
        self.cluster_tools.append(self.usearch)

    def run_orthofinder(self):
        """
        用orthofinder进行聚类
        :return:
        """
        self.logger.info("正在用orthofinder聚类软件进行聚类")
        self.orthofinder = self.add_tool("bac_comp_genome.pan_orthofinder")
        all_dir = self.convert_format.output_dir
        query_fasta = self.get_faa(all_dir)
        self.orthofinder.set_options({
            "infile_dir": query_fasta,
            "identity": self.option("identity"),
            "coverage": self.option("coverage"),
            "inflation": self.option("inflation"),
            "species": False,
            })
        self.orthofinder.on('end', self.set_output)
        self.orthofinder.run()
        self.cluster_tools.append(self.orthofinder)

    def run_orthomcl(self):
        """
        用orthomcl进行聚类
        :return:
        """
        self.logger.info("正在用orthomcl聚类软件进行聚类")
        self.orthomcl = self.add_tool("bac_comp_genome.pan_orthomcl")
        #compare = self.work_dir + '/cluster_tmp'
        coverage = self.option("coverage")*100
        identity = self.option("identity")*100
        all_dir = self.convert_format.output_dir
        query_fasta = self.get_faa(all_dir)
        self.orthomcl.set_options({
            "fasta_dir": query_fasta,
            "pi_cutoff": identity,
            "pmatch_cutoff": coverage,
            "inflation": self.option("inflation"),
            })
        self.orthomcl.on('end', self.set_output)
        self.orthomcl.run()
        self.cluster_tools.append(self.orthomcl)

    def run_homologus(self):
        """
        用Get_homologus进行聚类
        :return:
        """
        self.logger.info("正在用GET_HOMOLOGUS聚类软件进行聚类")
        self.homologus = self.add_tool("bac_comp_genome.pan_homologus")
        #compare = self.work_dir + '/cluster_tmp'
        all_dir = self.convert_format.output_dir
        query_fasta = self.get_faa(all_dir)
        self.homologus.set_options({
            "infile_dir": query_fasta,
            "identity": self.option("identity"),
            "coverage": self.option("coverage"),
            "inflation": self.option("inflation"),
            "method": self.option("homologus_mehod"),
            })
        self.homologus.on('end', self.set_output)
        self.homologus.run()
        self.cluster_tools.append(self.homologus)

    def run_pgap(self):
        """
        用PGAP进行聚类
        :return:
        """
        self.logger.info("正在用PGAP聚类软件进行聚类")
        self.pgap = self.add_tool("bac_comp_genome.pan_pgap")
        #compare = self.work_dir + '/cluster_tmp'
        # query_fasta = self.convert_format.output_dir
        query_fasta = self.option("fasta_dir").prop['path']
        self.logger.info(self.convert_format.option('out_dir'))
        opts = ({
            "infile_dir": query_fasta,
            "identity": self.option("identity"),
            "coverage": self.option("coverage"),
            "method": self.option("pgap_mehod"),
            "inflation": self.option("inflation"),
            })
        if self.option("pgap_mehod") in ['GF']:
            opts["score"] = self.option("score")
            opts["evalue"] = self.option("evalue")
        elif self.option("pgap_mehod") in ['MP']:
            opts["local"] = self.option("local")
            opts["global"] = self.option("global")
        self.pgap.set_options(opts)
        self.pgap.on('end', self.set_output)
        self.pgap.run()
        self.cluster_tools.append(self.pgap)

    def run_roary(self):
        """
        用Roary进行聚类
        :return:
        """
        self.logger.info("正在用Roary聚类软件进行聚类")
        self.roary = self.add_tool("bac_comp_genome.pan_roary")
        #compare = self.work_dir + '/cluster_tmp'
        all_dir = self.convert_format.output_dir
        query_fasta = self.get_faa(all_dir)
        self.roary.set_options({
            "infile_dir": query_fasta,
            "identity": self.option("identity"),
            "inflation": self.option("inflation"),
            "coverage": self.option("coverage"),
            })
        self.roary.on('end', self.set_output)
        self.roary.run()
        self.cluster_tools.append(self.roary)

    def run_annotation(self):
        """
        根据clusterid和注释结果进行合并
        :return:
        """
        self.logger.info("开始对聚类和注释结果进行合并")
        cluster_path = self.get_cluster()
        self.pan_anno.set_options({
            "cluster": cluster_path,
            "anno_file": self.option("anno_file"),
            })
        self.pan_anno.on('end', self.set_output)
        self.pan_anno.run()

    def get_faa(self, all_dir):
        """
        从所有文件中获得faa文件夹
        :return:
        """
        faa_dir = os.path.join(self.work_dir, 'faa_dir')
        if os.path.exists(faa_dir):
            shutil.rmtree(faa_dir)
        os.mkdir(faa_dir)
        for file in os.listdir(all_dir):
            file_path = os.path.join(all_dir, file)
            new_file_path = os.path.join(faa_dir, file)
            if file.endswith('.faa'):
                os.link(file_path, new_file_path)
        return faa_dir

    def get_cluster(self):
        """
        获取聚类结果文件
        :return:
        """
        for tool in self.cluster_tools:
            for file in os.listdir(tool.output_dir):
                file_path = os.path.join(self.work_dir, file)
                if os.path.exists(file_path):
                    os.remove(file_path)
                os.link(os.path.join(tool.output_dir, file), file_path)
        cluster = os.path.join(self.work_dir,'homologues_cluster.xls')
        return cluster

    def set_output(self):
        """
        设置结果文件目录
        :return:
        """
        self.logger.info('设置结果目录')
        self.get_cluster()
        if os.path.exists(os.path.join(self.output_dir, 'homologues_cluster.xls')):
            os.remove(os.path.join(self.output_dir, 'homologues_cluster.xls'))
        os.link(os.path.join(self.work_dir, 'homologues_cluster.xls'), os.path.join(self.output_dir, 'homologues_cluster.xls'))
        self.option('out', os.path.join(self.output_dir, 'homologues_cluster.xls'))
        if os.path.exists(os.path.join(self.output_dir, 'cluster_anno.xls')):
            os.remove(os.path.join(self.output_dir, 'cluster_anno.xls'))
        #os.link(os.path.join(self.pan_anno.output_dir, 'cluster_anno.xls'), os.path.join(self.output_dir, 'cluster_anno.xls'))
        self.logger.info('设置结果目录成功')
        self.end()

    def run(self):
        """
        开始运行了
        :return:
        """
        self.logger.info("开始运行")
        super(PanClusterModule, self).run()
        if self.option("cluster_method") in ["cdhit"]:
            self.convert_format.on('end', self.run_cdhit)
            self.convert_data()
        elif self.option("cluster_method") in ["mmseq"]:
            self.convert_format.on('end', self.run_mmseq)
            self.convert_data()

        elif self.option("cluster_method") in ["usearch"]:
            self.convert_format.on('end', self.run_usearch)
            self.convert_data()

        elif self.option("cluster_method") in ['blast+mcl'] and self.option("method") in ['orthofinder']:
            self.convert_format.on('end', self.run_orthofinder)
            self.convert_data()

        elif self.option("cluster_method") in ['blast+mcl'] and self.option("method") in ['orthomcl']:
            self.convert_format.on('end', self.run_orthomcl)

        elif self.option("cluster_method") in ['blast+mcl'] and self.option("method") in ['get_homologus']:
            self.convert_format.on('end', self.run_homologus)
            self.convert_data()

        elif self.option("cluster_method") in ['blast+mcl'] and self.option("method") in ['pgap']:
            # self.convert_data()
            self.run_pgap()

        elif self.option("cluster_method") in ['cdhit+blast+mcl'] and self.option("method") in ['roary']:
            self.convert_format.on('end', self.run_roary)
            self.convert_data()

        else:
            self.set_error("没有输入正确的聚类方法")



    def end(self):
        """
        结束了
        :return:
        """
        super(PanClusterModule, self).end()
