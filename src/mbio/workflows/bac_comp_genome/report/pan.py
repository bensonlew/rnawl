# -*- coding: utf-8 -*-
#__author__ = 'qingchen.zhang'
import os,shutil
import re
import types
from biocluster.workflow import Workflow
from bson import ObjectId
from biocluster.core.exceptions import OptionError
import datetime
from biocluster.file import download, exists
from collections import defaultdict
import shutil
import gevent
from mbio.packages.bac_comp_genome.common_function import link_dir,link_file



class PanWorkflow(Workflow):
    def __init__(self, wsheet_object):
        """
        细菌比较基因组泛基因组分析的第一步：
        聚类、泛基因组大小计算和变异分析；
        第二步 对聚类结果进行分组方案的合并
        第三步 对聚类的结果进行venn图的分析
        :return:
        """
        self._sheet = wsheet_object
        super(PanWorkflow, self).__init__(wsheet_object)
        options = [
            {"name": "fasta_dir", "type": "string"},  # 输入fasta的dir文件
            {"name": "identity", "type": "float", "default": 0.5},  ##给出cdhit的参数identity
            {"name": "coverage", "type": "float", "default": 0.01},  # 给出cdhit的参数coverage
            {"name": "method", "type": "string", 'default':'pgap'}, #输入方法类型{'orthofinder', 'orthomcl', 'get_homologus', 'pgap', 'roary'}
            {"name": "cluster_method", "type": "string",'default': 'blast+mcl'}, #cdhit usearch mmseq blast+mcl
            {"name": "inflation", "type": "float", "default": 1.5},  # Markov Inflation Index
            {"name": "pangenome", "type": "string", "default": "homologus"}, ## 计算泛基因组大小的公式的方法
            {"name": "cal_method", "type": "string", "default": "core_Tettelin"}, ###计算方法
            {"name": "score", "type": "int", "default": 40},#Minimum score in blast
            {"name": "evalue", "type": "string", "default": "1e-10"},#Maximal E-value in blastall
            {"name": "local", "type": "float", "default": 0.25},#Minimum local alignment overlap in MP method
            {"name": "global", "type": "float", "default": 0.5},#Minimum global alignment overlap in MP method
            {"name": "group", "type": "infile", "format": "meta.otu.group_table"},#输入group表的id
            {"name": "group_detail", "type": "string"}, ## 用于导表group表
            {"name": "main_id", "type": "string"},#输入pan_venn表的id
            {"name": "update_info", "type": "string"},
            {"name": "pgap_mehod", "type": "string", "default": "GF"},  # GF for GeneFamily method,  and MP for MultiParanoid method
            {"name": "homologus_mehod", "type": "string", "default": "OMCL"},  # 三种方法BDBH、OMCL、OCOG
            {"name": "anno_file", "type": "infile", "format": "sequence.profile_table"},  # 物种注释总览表
            {"name": "params", "type": "string"},#兼容接口加的参数
        ]
        self.pangenome = self.add_tool("bac_comp_genome.pan_genome")
        self.panvariation = self.add_tool("bac_comp_genome.pan_variation")
        self.pan_anno = self.add_tool("bac_comp_genome.pan_anno")
        self.cluster = self.add_module("bac_comp_genome.pan_cluster")
        self.IMPORT_REPORT_DATA = True
        self.IMPORT_REPORT_AFTER_END = False
        self.all_tools = []
        self.add_option(options)
        self.set_options(self._sheet.options())
        self.sample = []

    def check_options(self):
        """
        进行参数二次检查
        :return:
        """
        self.logger.info("开始pan的参数检查")
        if not self.option("identity"):
            raise OptionError("必须设置输入的identity")
        if not self.option("fasta_dir"):
            raise OptionError("请提供输入序列的文件夹！")
        else:
            self.logger.info(self.option('fasta_dir'))
        return True

    def run_cluster(self):
        """
        聚类分析
        :return:
        """
        self.logger.info("开始用聚类软件进行聚类")
        #compare = self.work_dir + '/cluster_tmp'
        # file_list = os.listdir(self.option('fasta_dir').prop['path'])
        # query_fasta = self.option('fasta_dir').prop["path"]
        opts = ({
            "fasta_dir": self.fasta_dir,
            "identity": self.option("identity"),
            "coverage": self.option("coverage"),
            "method": self.option("method"),
            "cluster_method": self.option("cluster_method"),
            # "anno_file": self.option("anno_file"),
            })
        if self.option("inflation"):
            opts["inflation"] = self.option("inflation")
        if self.option("pgap_mehod") != "":
            opts['pgap_mehod'] = self.option("pgap_mehod")
        if self.option("homologus_mehod") != "":
            opts['homologus_mehod'] = self.option("homologus_mehod")
        if self.option("score") != 0:
            opts["score"] = self.option("score")
        if self.option("evalue") != "":
            opts["evalue"] = self.option("evalue")
        if self.option("local") != 0.0:
            opts["local"] = self.option("local")
        if self.option("global") != 0.0:
            opts["global"] = self.option("global")
        self.cluster.set_options(opts)
        self.cluster.run()

    def run_genome(self):
        """
        根据聚类结果进行计算泛基因组的大小
        :return:
        """
        self.logger.info("开始计算泛基因组大小")
        # query_fasta = self.option('fasta_dir').prop["path"]
        self.pangenome.set_options({
            "cluster": self.cluster.option("out"),
            "pangenome": self.option("pangenome"),
            "cal_method": self.option("cal_method"),
            "cal_newgene": "average",
            "infile_dir": self.fasta_dir,
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
            "infile_dir": self.fasta_dir,
            # "cal_method": self.option("cal_method"),
            })
        self.all_tools.append(self.panvariation)

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
        将计算大小、变异分析并行运行
        :return:
        """
        if len(os.listdir(self.fasta_dir)) > 2:## 防止样本数为1的情况
            self.run_genome()
        # self.run_variation()
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
        result_path = os.path.join(self.output_dir)
        if os.path.exists(result_path):
            shutil.rmtree(result_path)
        os.mkdir(result_path)
        cluster_dir = result_path
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
        if os.path.exists(os.path.join(self.cluster.output_dir, "homologues_cluster.xls")):
            os.link(os.path.join(self.cluster.output_dir, "homologues_cluster.xls"), os.path.join(cluster_dir,"homologues_cluster.xls"))
        if os.path.exists(os.path.join(self.cluster.output_dir, "cluster_anno.xls")):
            os.link(os.path.join(self.cluster.output_dir, "cluster_anno.xls"), os.path.join(cluster_dir,"cluster_anno.xls"))
        elif os.path.exists(os.path.join(self.pan_anno.output_dir, "cluster_anno.xls")):
            os.link(os.path.join(self.pan_anno.output_dir, "cluster_anno.xls"), os.path.join(cluster_dir,"cluster_anno.xls"))
        if os.path.exists(os.path.join(self.panvariation.output_dir, "CDS_variation.xls")):
            os.link(os.path.join(self.panvariation.output_dir, "CDS_variation.xls"), os.path.join(cluster_dir,"CDS_variation.xls"))
        if os.path.exists(os.path.join(self.panvariation.output_dir, "CDS_variation_analysis.xls")):
            os.link(os.path.join(self.panvariation.output_dir, "CDS_variation_analysis.xls"), os.path.join(cluster_dir,"CDS_variation_analysis.xls"))

        self.logger.info('设置结果目录成功')
        self.set_db()

    def set_db(self):
        """
        将数据导入mongo
        :return:
        """
        self.logger.info('正在写入mongo数据库')
        self.remote_dir = self._sheet.output
        cluster_dir = self.output_dir
        api_pan = self.api.api('bac_comp_genome.pan')

        cluster_data = os.path.join(cluster_dir, 'homologues_cluster.xls')
        # variation_data = os.path.join(cluster_dir, 'CDS_variation_analysis.xls')
        anno_file = os.path.join(cluster_dir, 'cluster_anno.xls')
        cluster_path = self.remote_dir + "homologues_cluster.xls"

        self.logger.info("正在导泛基因组分析的主表")
        pan_id = self.option("main_id")
        self.logger.info("正在导泛基因组分析的cluster结果表")
        api_pan.add_pan_detail(pan_id, cluster_data, anno_file, cluster_path=cluster_path)

        pan_graph = os.path.join(cluster_dir, 'pan_cluster_genome.xls')
        pan_detail = os.path.join(cluster_dir, 'pan_cluster.xls')
        core_graph = os.path.join(cluster_dir, 'core_cluster_genome.xls')
        core_detail = os.path.join(cluster_dir, 'core_cluster.xls')
        new_graph = os.path.join(cluster_dir, 'new_gene_cluster_genome.xls')
        new_detail = os.path.join(cluster_dir, 'new_gene_cluster.xls')

        self.logger.info("正在导公式的结果表")
        if os.path.exists(pan_graph) and os.path.exists(pan_detail):
            with open(pan_graph, 'r') as f:
                line_counts = len(f.readlines())
                if line_counts > 1:
                    api_pan.add_pan_genome(pan_id, 'pan', pan_graph, pan_detail)
        if os.path.exists(core_graph) and os.path.exists(core_detail):
            with open(core_graph, 'r') as f1:
                line_counts2= len(f1.readlines())
                if line_counts2 > 1:
                    api_pan.add_pan_genome(pan_id, 'core', core_graph, core_detail)
        if os.path.exists(new_graph) and os.path.exists(new_detail):
            with open(new_graph, 'r') as f2:
                line_counts3= len(f2.readlines())
                if line_counts3 > 1:
                    api_pan.add_pan_genome(pan_id, 'newgene', new_graph, new_detail)

        pgap_formula = os.path.join(cluster_dir, 'pgap_formula.xls')
        if os.path.exists(pgap_formula):
            api_pan.add_pan_formula(pan_id, pgap_formula)

        self.logger.info("正在导变异分析的统计结果表")
        variation_path = os.path.join(cluster_dir, 'CDS_variation_analysis.xls')
        var_number = 0
        if os.path.exists(variation_path):
            with open(variation_path, 'r') as var:
                for line in var:
                    var_number += 1
            if var_number > 1:
                api_pan.add_pan_variation_stat(pan_id, variation_path)
        self.end()

    def get_fasta_dir(self):
        """
        根据group和输入的文件夹获取此次分析的样本夹
        :return:
        """
        self.fasta_dir = os.path.join(self.work_dir, "fasta_dir")
        if os.path.exists(self.fasta_dir):
            shutil.rmtree(self.fasta_dir)
        os.mkdir(self.fasta_dir)
        with open(self.option("group").prop["path"], "r") as f:
            f.readline()
            for line in f:
                line = line.strip().split("\t")
                sample_name = line[0]
                for i in ["_CDS.fna", "_CDS.faa"]:
                    if os.path.exists(self.option("fasta_dir") + "/" + sample_name + "/" + sample_name + i):
                        self.logger.info("fasta_dir: {}".format(self.option("fasta_dir")))
                        if "_CDS.fnn" == i:
                            download(self.option("fasta_dir") + "/" + sample_name + "/" + sample_name + i, self.fasta_dir+ "/" + sample_name + "_CDS.fna")
                        else:
                            download(self.option("fasta_dir") + "/" + sample_name + "/" + sample_name + i, self.fasta_dir + "/" + sample_name + i)
                    else:
                        download(self.option("fasta_dir") + "/" + sample_name + "/" + sample_name + i, self.fasta_dir+ "/" + sample_name + i)

    def run(self):
        """
        开始运行了
        :return:
        """
        self.logger.info("开始运行")
        self.get_fasta_dir()
        self.cluster.on('end', self.run_all_tool)
        self.run_cluster()
        super(PanWorkflow, self).run()

    def end(self):
        """
        结束了
        :return:
        """
        self.logger.info("开始结果文件上传")
        repaths = [
            [".", "", "泛基因组分析结果目录",0,""],
            ["homologues_cluster.xls", "xls", "同源基因结果表",0,""],
            ["CDS_variation_analysis.xls", "xls", "同源基因变异分析结果表",0,""],
            ["CDS_variation.xls", "xls", "同源基因变异分析统计表",0,""],
            ["cluster_anno.xls", "xls", "同源基因注释结果表",0,""],
            ["pan_cluster_genome.xls", "xls", "pangenomes计算公式结果表",0,""],
            ["pan_cluster.xls", "xls", "pangenomes曲线结果表",0,""],
            ["core_cluster_genome.xls", "xls", "coregene计算公式结果表",0,""],
            ["core_cluster.xls", "xls", "coregene曲线结果表",0,""],
            ["new_gene_cluster_genome.xls", "xls", "newgene计算公式结果表",0,""],
            ["new_gene_cluster.xls", "xls", "newgene曲线结果表",0,""],
            ["pgap_formula.xls", "xls", "pangenomes计算公式(pgap)结果表",0,""],

        ]
        sdir = self.add_upload_dir(self.output_dir)
        sdir.add_relpath_rules(repaths)
        super(PanWorkflow, self).end()
