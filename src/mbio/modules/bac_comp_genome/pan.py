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


class PanModule(Module):
    def __init__(self, work_id):
        """
        细菌比较基因组泛基因组分析的第一步：
        聚类、泛基因组大小计算和变异分析；
        第二步 对聚类结果进行分组方案的合并
        第三步 对聚类的结果进行venn图的分析
        :return:
        """
        super(PanModule, self).__init__(work_id)
        options = [
            {"name": "fasta_dir", "type": "infile", "format": "bac_comp_genome.input_dir"},  # 输入fasta的dir文件包含核酸和蛋白
            {"name": "identity", "type": "float", "default": 0.5},  ##给出cdhit的参数identity
            {"name": "coverage", "type": "float", "default": 0.0},  # 给出cdhit的参数coverage
            {"name": "method", "type": "string"}, #输入方法类型{'orthofinder', 'orthomcl', 'get_homologus', 'pgap', 'roary'}
            {"name": "cluster_method", "type": "string"}, #cdhit usearch mmseq blast+mcl
            {"name": "inflation", "type": "float", "default": 1.5},  # Markov Inflation Index
            {"name": "pangenome", "type": "string", "default": "homologus"}, ## 计算泛基因组大小的公式的方法
            {"name": "cal_method", "type": "string", "default": "core_Tettelin"}, ###计算方法
            {"name": "category", "type": "string"}, #输入分类方案类型{'core'：{"min":1, "max": 1}, 'dispensable'：{"min":0.15, "max": 0.95}, 'soft_core'{"min":0.95, "max": 1}, 'unique'{"min":0, "max": 0.15}}min能取到，max取不到
            {"name": "percent", "type":"bool", "default": False}, #false 表示按个数，true表示按百分比
            {"name": "category_name", "type":"string"}, #输入分组方案的名称
            {"name": "group_table", "type": "infile", "format": "meta.otu.group_table"},#输入group表选择哪些样本进行导表和计算
            {"name": "anno_file", "type": "infile", "format": "sequence.profile_table"},  # 物种注释总览表
        ]
        self.pangenome = self.add_tool("bac_comp_genome.pan_genome")
        self.panvariation = self.add_tool("bac_comp_genome.pan_variation")
        self.category = self.add_tool("bac_comp_genome.pan_category")
        self.venn = self.add_tool("bac_comp_genome.pan_venn")
        self.cluster = self.add_module("bac_comp_genome.pan_cluster")
        self.pan_anno = self.add_tool("bac_comp_genome.pan_anno")
        self.all_tools = []
        self.add_option(options)

    def check_options(self):
        """
        进行参数二次检查
        :return:
        """
        if not self.option("identity"):
            raise OptionError("必须设置输入的identity")
        return True

    def run_cluster(self):
        """
        聚类分析
        :return:
        """
        self.logger.info("开始用聚类软件进行聚类")
        #compare = self.work_dir + '/cluster_tmp'
        file_list = os.listdir(self.option('fasta_dir').prop['path'])
        file_num = len(file_list) / 2
        if file_num >= 60:
            cluster_method = "cdhit+blast+mcl"
            method = "roary"
        else:
            cluster_method = "blast+mcl"
            method = "pgap"
        query_fasta = self.option('fasta_dir').prop["path"]
        self.cluster.set_options({
            "fasta_dir": query_fasta,
            "identity": self.option("identity"),
            "coverage": self.option("coverage"),
            "inflation": self.option("inflation"),
            "method": method,
            "cluster_method": cluster_method,
            # "anno_file": self.option("anno_file"),
            })
        self.cluster.run()

    def run_genome(self):
        """
        根据聚类结果进行计算泛基因组的大小
        :return:
        """
        self.logger.info("开始计算泛基因组大小")
        query_fasta = self.option('fasta_dir').prop["path"]
        self.pangenome.set_options({
            "cluster": self.cluster.option("out"),
            "pangenome": self.option("pangenome"),
            "cal_method": self.option("cal_method"),
            "cal_newgene": "average",
            "infile_dir": query_fasta,
            })
        self.all_tools.append(self.pangenome)

    def run_variation(self):
        """
        对聚类结果进行变异分析
        :return:
        """
        self.logger.info("开始进行变异分析")
        self.panvariation.set_options({
            "cluster": self.cluster.option("out"),
            "infile_dir": self.option("fasta_dir"),
            # "cal_method": self.option("cal_method"),
            })
        self.all_tools.append(self.panvariation)

    def run_category(self):
        """
        根据分类方案计算整理
        :return:
        """
        self.logger.info("开始进行分组方案合并")
        opts = ({
            "cluster": self.cluster.option("out"),
            "category": self.option("category"),
            "percent": self.option("percent"),
            "category_name": self.option("category_name"),
            })
        # if self.option("group_table").is_set:
        #     opts["group_table"] = self.option("group_table")
        self.category.set_options(opts)
        self.all_tools.append(self.category)

    def run_venn(self):
        """
        Venn图的分析
        :return:
        """
        self.logger.info("开始进行venn图分析")
        self.venn.set_options({
            "cluster": self.cluster.option("out"),
            "group_table": self.option("group_table"),
            "version": "v1"
            })
        self.all_tools.append(self.venn)

    def run_annotation(self):
        """
        根据clusterid和注释结果进行合并
        :return:
        """
        self.logger.info("开始对聚类和注释结果进行合并")
        # cluster_path = self.get_cluster()
        self.pan_anno.set_options({
            "cluster": self.cluster.option("out"),
            "anno_file": self.option("anno_file"),
            })
        # self.pan_anno.on('end', self.set_output)
        self.all_tools.append(self.pan_anno)

    def run_all_tool(self):
        """
        将计算大小、变异分析、分组方案和venn图分析并行运行
        :return:
        """
        self.run_genome()
        # self.run_variation()
        self.run_category()
        self.run_venn()
        self.run_annotation()
        self.on_rely(self.all_tools, self.set_output)
        for tool in self.all_tools:
            tool.run()
            gevent.sleep(0)

    def set_output(self):
        """
        设置结果文件目录
        :return:
        """
        self.logger.info('设置结果目录')
        result_path = os.path.join(self.output_dir, 'pangenomes')
        if os.path.exists(result_path):
            shutil.rmtree(result_path)
        os.mkdir(result_path)
        cluster_dir = os.path.join(result_path, "cluster")
        category_dir = os.path.join(result_path, "category")
        venn_dir = os.path.join(result_path, "venn")
        if os.path.exists(cluster_dir):
            shutil.rmtree(cluster_dir)
        if os.path.exists(category_dir):
            shutil.rmtree(category_dir)
        if os.path.exists(venn_dir):
            shutil.rmtree(venn_dir)
        os.mkdir(cluster_dir)
        os.mkdir(category_dir)
        os.mkdir(venn_dir)
        # link_dir(self.pangenome.output_dir, cluster_dir)
        for file in os.listdir(self.pangenome.output_dir):
            if file.endswith('.log'):
                if file.startswith('core'):
                    os.link(os.path.join(self.pangenome.output_dir, file), cluster_dir + "/core_cluster_genome.xls")
                elif file.startswith('pan'):
                    os.link(os.path.join(self.pangenome.output_dir, file), cluster_dir + "/pan_cluster_genome.xls")
                else:
                    os.link(os.path.join(self.pangenome.output_dir, file), cluster_dir +
                            "/new_gene_cluster_genome.xls")
            else:
                os.link(os.path.join(self.pangenome.output_dir, file), cluster_dir+"/" + file)
        os.link(os.path.join(self.cluster.output_dir, "homologues_cluster.xls"), os.path.join(cluster_dir,"homologues_cluster.xls"))
        if os.path.exists(os.path.join(self.cluster.output_dir, "cluster_anno.xls")):
            os.link(os.path.join(self.cluster.output_dir, "cluster_anno.xls"), os.path.join(cluster_dir,"cluster_anno.xls"))
        elif os.path.exists(os.path.join(self.pan_anno.output_dir, "cluster_anno.xls")):
            os.link(os.path.join(self.pan_anno.output_dir, "cluster_anno.xls"), os.path.join(cluster_dir,"cluster_anno.xls"))
        # os.link(os.path.join(self.panvariation.output_dir, "CDS_variation.xls"), os.path.join(cluster_dir,"CDS_variation.xls"))
        # os.link(os.path.join(self.panvariation.output_dir, "CDS_variation_analysis.xls"), os.path.join(cluster_dir,"CDS_variation_analysis.xls"))

        link_dir(self.category.output_dir,category_dir)
        link_dir(self.venn.output_dir, venn_dir)
        # os.link(os.path.join(self.cluster.output_dir, "merge_category.xls"), category_dir)
        # os.link(os.path.join(self.cluster.output_dir, "venn_stat.xls"), venn_dir)
        self.logger.info('设置结果目录成功')
        self.end()

    def run(self):
        """
        开始运行了
        :return:
        """
        super(PanModule, self).run()
        self.logger.info("开始运行")
        self.cluster.on('end', self.run_all_tool)
        self.run_cluster()



    def end(self):
        """
        结束了
        :return:
        """
        super(PanModule, self).end()
