# -*- coding: utf-8 -*-

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
            {"name": "coverage", "type": "float", "default": 0.0000001},  # 给出cdhit的参数coverage
            {"name": "method", "type": "string"}, #输入方法类型{'orthofinder', 'orthomcl', 'get_homologus', 'pgap', 'roary'}
            #{"name": "cluster_method", "type": "string"}, #cdhit usearch mmseq blast+mcl
            {"name": "inflation", "type": "float", "default": 1.5},  # Markov Inflation Index
            {"name": "pangenome", "type": "string", "default": "homologus"}, ## 计算泛基因组大小的公式的方法
            {"name": "cal_method", "type": "string", "default": "core_Tettelin"}, ###计算方法
            {"name": "homologus_mehod", "type": "string", "default": "OMCL"},  # 三种方法BDBH、OMCL、OCOG
            {"name": "pgap_mehod", "type": "string", "default": "GF"},  # GF MP
            {"name": "score", "type": "int", "default": 40},
            {"name": "evalue", "type": "string", "default": '1e-10'},
            {"name": "local", "type": "float", "default": 0.25},#Minimum local alignment overlap in MP method
            {"name": "global", "type": "float", "default": 0.5},
            {"name": "group_table", "type": "infile", "format": "meta.otu.group_table"},#输入group表选择哪些样本进行导表和计算 例如： #sample group
            {"name": "cal_pgap", "type": "string", "default": "false"},
        ]
        self.pangenome = self.add_tool("tool_lab.pan_genome")
        self.venn = self.add_tool("tool_lab.pan_venn")
        self.cluster = self.add_module("tool_lab.pan_cluster")
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
        cluster_method = "blast+mcl"
        query_fasta = self.option('fasta_dir').prop["path"]
        method = self.option('method')

        opts = {
            "fasta_dir": query_fasta,
            "method": method,
            "cluster_method": cluster_method
        }
        if method == 'get_homologus':
            opts['homologus_mehod'] = self.option("homologus_mehod")
            opts["identity"] = self.option("identity")
            opts["coverage"] =  self.option("coverage")
            opts["inflation"] = self.option("inflation")

        elif method == 'pgap':
            opts["pgap_mehod"] = self.option("pgap_mehod")
            opts["identity"] = self.option("identity")
            opts["coverage"] =  self.option("coverage")
            opts["inflation"] = self.option("inflation")

            if self.option("pgap_mehod") == 'GF':
                opts["score"] = self.option("score")
                opts["evalue"] = self.option("evalue")
            elif self.option("pgap_mehod") == 'MP':
                opts["local"] = self.option("local")
                opts["global"] = self.option("global")

        self.cluster.set_options(opts)
        self.cluster.run()

    def run_genome(self):
        """
        根据聚类结果进行计算泛基因组的大小
        :return:
        """
        self.logger.info("开始计算泛基因组大小")
        query_fasta = self.option('fasta_dir').prop["path"]
        opt = {
            "cluster": self.cluster.option("out"),
            "pangenome": self.option("pangenome"),
            "cal_method": self.option("cal_method"),
            "cal_newgene": "average",
            "infile_dir": query_fasta
            }
        if self.option('cal_pgap') == 'true':
            opt['cal_method'] = 'core_both'
        self.pangenome.set_options(opt)
        self.all_tools.append(self.pangenome)

    def run_venn(self):
        """
        根据同源基因Venn图的分析
        :return:
        """
        self.logger.info("开始进行venn图分析")
        self.venn.set_options({
            "cluster": self.cluster.option("out"),
            "group_table": self.option("group_table"),
            "version": "v1"
            })
        self.all_tools.append(self.venn)


    def run_all_tool(self):
        """
        将计算大小、变异分析、分组方案和cluster的venn图分析并行运行
        :return:
        """
        self.run_genome()
        self.run_venn()
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
        venn_dir = os.path.join(result_path, "venn")
        if os.path.exists(cluster_dir):
            shutil.rmtree(cluster_dir)
        if os.path.exists(venn_dir):
            shutil.rmtree(venn_dir)
        os.mkdir(cluster_dir)
        os.mkdir(venn_dir)
        for file in os.listdir(self.pangenome.output_dir):
            if file.endswith('.log'):
                if file.startswith('core'):
                    if os.path.exists(cluster_dir + "/core_cluster_genome.xls"):
                        os.remove(cluster_dir + "/core_cluster_genome.xls")
                    os.link(os.path.join(self.pangenome.output_dir, file), cluster_dir + "/core_cluster_genome.xls")
                elif file.startswith('pan'):
                    if os.path.exists(cluster_dir + "/pan_cluster_genome.xls"):
                        os.remove(cluster_dir + "/pan_cluster_genome.xls")
                    os.link(os.path.join(self.pangenome.output_dir, file), cluster_dir + "/pan_cluster_genome.xls")
                else:
                    if os.path.exists(cluster_dir + "/new_gene_cluster_genome.xls"):
                        os.remove(cluster_dir + "/new_gene_cluster_genome.xls")

                    os.link(os.path.join(self.pangenome.output_dir, file), cluster_dir +
                            "/new_gene_cluster_genome.xls")
            else:
                os.link(os.path.join(self.pangenome.output_dir, file), cluster_dir + "/" + file)
        if os.path.exists(os.path.join(cluster_dir, "homologues_cluster.xls")):
            os.remove(os.path.join(cluster_dir, "homologues_cluster.xls"))

        os.link(os.path.join(self.cluster.output_dir, "homologues_cluster.xls"), os.path.join(cluster_dir, "homologues_cluster.xls"))
        link_dir(self.venn.output_dir, venn_dir)
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
